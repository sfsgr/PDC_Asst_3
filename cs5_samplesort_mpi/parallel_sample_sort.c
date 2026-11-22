#include "common.h"
#include <stdio.h>
#include <string.h>
#include <mpi.h>

static uint32_t*lower_bound(uint32_t*begin,uint32_t*end,uint32_t value) {
  while(begin<end) {
    uint32_t*middle=begin+(end-begin)/2;
    if(*middle<=value) {
      begin=middle+1;
    } else {
      end=middle;
    }
  }
  return begin;
}

static int verify_global_sorted(uint32_t*data,size_t count,int rank,int ranks,MPI_Comm comm) {
  int local_ok=verify_local_sorted(data,count);
  uint32_t local_min=count?data[0]:UINT32_MAX;
  uint32_t local_max=count?data[count-1]:0;
  uint32_t previous_max=0;
  if(rank>0) {
    MPI_Recv(&previous_max,1,MPI_UINT32_T,rank-1,100,comm,MPI_STATUS_IGNORE);
  }
  if(rank<ranks-1) {
    MPI_Send(&local_max,1,MPI_UINT32_T,rank+1,100,comm);
  }
  if(rank>0&&count&&previous_max>local_min) {
    local_ok=0;
  }
  int global_ok=0;
  MPI_Allreduce(&local_ok,&global_ok,1,MPI_INT,MPI_MIN,comm);
  return global_ok;
}

int main(int argc,char**argv) {
  MPI_Init(&argc,&argv);
  int rank,ranks;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  size_t n=(argc>1)?strtoull(argv[1],NULL,10):10000000;
  int distribution=(argc>2)?atoi(argv[2]):0;
  int oversample=(argc>3)?atoi(argv[3]):1;
  if(oversample<1) {
    oversample=1;
  }
  size_t base=n/ranks;
  size_t remainder=n%ranks;
  size_t local_n=base+(rank<(int)remainder);
  uint32_t*local=(uint32_t*)malloc(local_n*sizeof(uint32_t));
  if(!local) {
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  if(distribution==0) {
    generate_uniform_data(local,local_n,42ULL+(uint64_t)rank*10007ULL);
  } else {
    generate_skewed_data(local,local_n,42ULL+(uint64_t)rank*10007ULL);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  double start=MPI_Wtime();
  double phase1_start=MPI_Wtime();
  qsort(local,local_n,sizeof(uint32_t),uint32_compare);
  double phase1=MPI_Wtime()-phase1_start;

  double phase2_start=MPI_Wtime();
  int sample_count=(ranks-1)*oversample;
  if(sample_count<1) {
    sample_count=1;
  }
  uint32_t*samples=(uint32_t*)malloc((size_t)sample_count*sizeof(uint32_t));
  for(int i=0;i<sample_count;i++) {
    size_t index=(size_t)(i+1)*local_n/(sample_count+1);
    if(index>=local_n) {
      index=local_n-1;
    }
    samples[i]=local[index];
  }
  uint32_t*gathered=rank==0?(uint32_t*)malloc((size_t)sample_count*ranks*sizeof(uint32_t)):NULL;
  MPI_Gather(samples,sample_count,MPI_UINT32_T,gathered,sample_count,MPI_UINT32_T,0,MPI_COMM_WORLD);
  uint32_t*splitters=(uint32_t*)malloc((size_t)(ranks-1)*sizeof(uint32_t));
  if(rank==0) {
    size_t total=(size_t)sample_count*ranks;
    qsort(gathered,total,sizeof(uint32_t),uint32_compare);
    for(int i=0;i<ranks-1;i++) {
      size_t index=(size_t)(i+1)*total/ranks;
      if(index>=total) {
        index=total-1;
      }
      splitters[i]=gathered[index];
    }
  }
  MPI_Bcast(splitters,ranks-1,MPI_UINT32_T,0,MPI_COMM_WORLD);
  double phase2=MPI_Wtime()-phase2_start;

  double phase3_start=MPI_Wtime();
  int*send_counts=(int*)calloc(ranks,sizeof(int));
  int*receive_counts=(int*)calloc(ranks,sizeof(int));
  int*send_displacements=(int*)malloc(ranks*sizeof(int));
  int*receive_displacements=(int*)malloc(ranks*sizeof(int));
  size_t previous=0;
  for(int i=0;i<ranks-1;i++) {
    uint32_t*position=lower_bound(local+previous,local+local_n,splitters[i]);
    size_t index=(size_t)(position-local);
    send_counts[i]=(int)(index-previous);
    previous=index;
  }
  send_counts[ranks-1]=(int)(local_n-previous);
  double phase3=MPI_Wtime()-phase3_start;

  double phase4_start=MPI_Wtime();
  MPI_Alltoall(send_counts,1,MPI_INT,receive_counts,1,MPI_INT,MPI_COMM_WORLD);
  int total_receive=0;
  for(int i=0;i<ranks;i++) {
    send_displacements[i]=(i==0)?0:send_displacements[i-1]+send_counts[i-1];
    receive_displacements[i]=(i==0)?0:receive_displacements[i-1]+receive_counts[i-1];
    total_receive+=receive_counts[i];
  }
  double phase4=MPI_Wtime()-phase4_start;

  double phase5_start=MPI_Wtime();
  uint32_t*received=(uint32_t*)malloc((size_t)total_receive*sizeof(uint32_t));
  MPI_Alltoallv(local,send_counts,send_displacements,MPI_UINT32_T,
                received,receive_counts,receive_displacements,MPI_UINT32_T,MPI_COMM_WORLD);
  double phase5=MPI_Wtime()-phase5_start;

  double phase6_start=MPI_Wtime();
  qsort(received,(size_t)total_receive,sizeof(uint32_t),uint32_compare);
  double phase6=MPI_Wtime()-phase6_start;
  MPI_Barrier(MPI_COMM_WORLD);
  double total=MPI_Wtime()-start;

  int sorted=verify_global_sorted(received,(size_t)total_receive,rank,ranks,MPI_COMM_WORLD);
  unsigned long long local_count=(unsigned long long)total_receive;
  unsigned long long global_count=0;
  MPI_Allreduce(&local_count,&global_count,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
  int count_ok=global_count==(unsigned long long)n;

  double phases[7]={total,phase1,phase2,phase3,phase4,phase5,phase6};
  double maximum[7];
  MPI_Reduce(phases,maximum,7,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  int*counts=(int*)malloc(ranks*sizeof(int));
  int final_count=total_receive;
  MPI_Gather(&final_count,1,MPI_INT,counts,1,MPI_INT,0,MPI_COMM_WORLD);

  if(rank==0) {
    int minimum=counts[0],maximum_count=counts[0];
    for(int i=1;i<ranks;i++) {
      if(counts[i]<minimum) {
        minimum=counts[i];
      }
      if(counts[i]>maximum_count) {
        maximum_count=counts[i];
      }
    }
    printf("======================================================\n");
    printf("Case Study 5: Distributed Sample Sort\n");
    printf("MPI Processes: %d\n",ranks);
    printf("Total Keys: %zu\n",n);
    printf("Distribution: %s\n",distribution==0?"Uniform":"Skewed");
    printf("Total Time: %.4f s\n",maximum[0]);
    printf("Throughput: %.2f Million keys/s\n",(n/1e6)/maximum[0]);
    printf("Local Sort: %.4f s\n",maximum[1]);
    printf("Sampling: %.4f s\n",maximum[2]);
    printf("Partition: %.4f s\n",maximum[3]);
    printf("Alltoall Counts: %.4f s\n",maximum[4]);
    printf("Alltoallv Data: %.4f s\n",maximum[5]);
    printf("Final Sort: %.4f s\n",maximum[6]);
    printf("Min Elements: %d\n",minimum);
    printf("Max Elements: %d\n",maximum_count);
    printf("Imbalance: %.3f\n",(double)maximum_count/((double)n/ranks));
    printf("Verification: Sorted=%s Count=%s\n",sorted?"PASSED":"FAILED",count_ok?"PASSED":"FAILED");
  }
  free(counts);
  free(received);
  free(receive_displacements);
  free(send_displacements);
  free(receive_counts);
  free(send_counts);
  free(splitters);
  if(gathered) {
    free(gathered);
  }
  free(samples);
  free(local);
  MPI_Finalize();
  return sorted&&count_ok?0:1;
}
