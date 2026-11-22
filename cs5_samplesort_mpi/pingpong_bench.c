#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc,char**argv) {
  MPI_Init(&argc,&argv);
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if(size<2) {
    if(rank==0) {
      fprintf(stderr,"ping-pong requires at least 2 ranks\n");
    }
    MPI_Finalize();
    return 1;
  }
  size_t sizes[]={8,64,512,4096,32768,131072,524288,1048576,4194304,16777216};
  int iterations=100;
  int warmup=10;
  size_t size_count=sizeof(sizes)/sizeof(sizes[0]);
  double*measured=(double*)malloc(size_count*sizeof(double));
  double*message_sizes=(double*)malloc(size_count*sizeof(double));
  if(!measured||!message_sizes) {
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  if(rank==0) {
    printf("======================================================\n");
    printf("Case Study 5: MPI Ping-Pong Benchmark\n");
    printf("Ranks: 0 <-> 1\n");
  }
  for(size_t s=0;s<sizeof(sizes)/sizeof(sizes[0]);s++) {
    size_t bytes=sizes[s];
    char*buffer=(char*)malloc(bytes);
    if(!buffer) {
      MPI_Abort(MPI_COMM_WORLD,1);
    }
    for(int i=0;i<warmup;i++) {
      if(rank==0) {
        MPI_Send(buffer,bytes,MPI_BYTE,1,0,MPI_COMM_WORLD);
        MPI_Recv(buffer,bytes,MPI_BYTE,1,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      } else if(rank==1) {
        MPI_Recv(buffer,bytes,MPI_BYTE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Send(buffer,bytes,MPI_BYTE,0,1,MPI_COMM_WORLD);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double start=MPI_Wtime();
    for(int i=0;i<iterations;i++) {
      if(rank==0) {
        MPI_Send(buffer,bytes,MPI_BYTE,1,0,MPI_COMM_WORLD);
        MPI_Recv(buffer,bytes,MPI_BYTE,1,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      } else if(rank==1) {
        MPI_Recv(buffer,bytes,MPI_BYTE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Send(buffer,bytes,MPI_BYTE,0,1,MPI_COMM_WORLD);
      }
    }
    double elapsed=MPI_Wtime()-start;
    if(rank==0) {
      double one_way=elapsed/(2.0*iterations);
      double bandwidth=(bytes/(1024.0*1024.0))/one_way;
      printf("Size=%zu bytes OneWay=%.3f us Bandwidth=%.2f MB/s\n",bytes,one_way*1e6,bandwidth);
      measured[s]=one_way;
      message_sizes[s]=(double)bytes;
    }
    free(buffer);
  }
  if(rank==0) {
    double sum_x=0.0;
    double sum_y=0.0;
    double sum_xx=0.0;
    double sum_xy=0.0;
    for(size_t i=0;i<size_count;i++) {
      sum_x+=message_sizes[i];
      sum_y+=measured[i];
      sum_xx+=message_sizes[i]*message_sizes[i];
      sum_xy+=message_sizes[i]*measured[i];
    }
    double count=(double)size_count;
    double beta=(count*sum_xy-sum_x*sum_y)/(count*sum_xx-sum_x*sum_x);
    double alpha=(sum_y-beta*sum_x)/count;
    if(alpha<0.0) {
      alpha=measured[0];
    }
    printf("Hockney Alpha: %.3f us\n",alpha*1e6);
    printf("Hockney Beta: %.4f ns/byte\n",beta*1e9);
    printf("Hockney Bandwidth: %.2f GB/s\n",(1.0/beta)/(1024.0*1024.0*1024.0));
  }
  free(measured);
  free(message_sizes);
  MPI_Finalize();
  return 0;
}
