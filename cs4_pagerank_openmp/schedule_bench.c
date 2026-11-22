#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc,char**argv) {
  int scale=(argc>1)?atoi(argv[1]):16;
  size_t target=(argc>2)?strtoull(argv[2],NULL,10):1500000;
  int threads=(argc>3)?atoi(argv[3]):8;
  omp_set_dynamic(0);
  omp_set_num_threads(threads);
  CSRGraph graph=generate_rmat_graph(scale,target,42ULL);
  double*rank=(double*)malloc((size_t)graph.n*sizeof(double));
  double*next=(double*)malloc((size_t)graph.n*sizeof(double));
  if(!rank||!next) {
    fprintf(stderr,"allocation failed\n");
    return 1;
  }
  for(int i=0;i<graph.n;i++) {
    rank[i]=1.0/graph.n;
  }
  printf("======================================================\n");
  printf("Case Study 4: OpenMP Schedule Benchmark\n");
  printf("Vertices: %d Edges: %zu Threads: %d\n",graph.n,graph.edge_count,threads);
  printf("In-Degree Distribution:\n");
  int degree_limits[]={0,1,2,4,8,16,32,64,128,256,512,1024};
  int degree_upper[]={0,1,3,7,15,31,63,127,255,511,1023,graph.n};
  for(int bin=0;bin<12;bin++) {
    int vertex_count=0;
    size_t edge_count=0;
    for(int v=0;v<graph.n;v++) {
      int degree=graph.in_degree[v];
      if(degree>=degree_limits[bin]&&degree<=degree_upper[bin]) {
        vertex_count++;
        edge_count+=(size_t)degree;
      }
    }
    printf("%d-%d: vertices=%d edges=%zu\n",degree_limits[bin],degree_upper[bin],vertex_count,edge_count);
  }
  double*thread_time=(double*)calloc((size_t)threads,sizeof(double));
  size_t*thread_edges=(size_t*)calloc((size_t)threads,sizeof(size_t));
  for(int iteration=0;iteration<10;iteration++) {
    #pragma omp parallel
    {
      int tid=omp_get_thread_num();
      size_t local_edges=0;
      double start_thread=omp_get_wtime();
      #pragma omp for schedule(static)
      for(int v=0;v<graph.n;v++) {
        double incoming=0.0;
        for(size_t e=graph.in_row_ptr[v];e<graph.in_row_ptr[v+1];e++) {
          int u=graph.in_col_idx[e];
          incoming+=rank[u]/graph.out_degree[u];
          local_edges++;
        }
        next[v]=0.15/graph.n+0.85*incoming;
      }
      thread_time[tid]+=omp_get_wtime()-start_thread;
      thread_edges[tid]=local_edges;
    }
  }
  double max_time=thread_time[0];
  double mean_time=0.0;
  size_t min_edges=thread_edges[0];
  size_t max_edges=thread_edges[0];
  for(int t=0;t<threads;t++) {
    if(thread_time[t]>max_time) {
      max_time=thread_time[t];
    }
    mean_time+=thread_time[t];
    if(thread_edges[t]<min_edges) {
      min_edges=thread_edges[t];
    }
    if(thread_edges[t]>max_edges) {
      max_edges=thread_edges[t];
    }
  }
  mean_time/=threads;
  printf("Static imbalance ratio: %.3f (slowest/mean)\n",max_time/mean_time);
  printf("Static work imbalance: %.2f (max edges/min edges)\n",(double)max_edges/(min_edges?min_edges:1));
  for(int t=0;t<threads;t++) {
    printf("Thread %d: time=%.6f edges=%zu\n",t,thread_time[t],thread_edges[t]);
  }
  printf("Schedule Chunk Time(s) Throughput(Medges/s)\n");
  const char*names[]={ "static","static","dynamic","dynamic","guided","guided" };
  int types[]={0,0,1,1,2,2};
  int chunks[]={64,256,64,256,64,256};
  for(int config=0;config<6;config++) {
    double start=omp_get_wtime();
    for(int iteration=0;iteration<15;iteration++) {
      if(types[config]==0) {
        #pragma omp parallel for schedule(static,chunks[config])
        for(int v=0;v<graph.n;v++) {
          double incoming=0.0;
          for(size_t e=graph.in_row_ptr[v];e<graph.in_row_ptr[v+1];e++) {
            int u=graph.in_col_idx[e];
            incoming+=rank[u]/graph.out_degree[u];
          }
          next[v]=0.15/graph.n+0.85*incoming;
        }
      } else if(types[config]==1) {
        #pragma omp parallel for schedule(dynamic,chunks[config])
        for(int v=0;v<graph.n;v++) {
          double incoming=0.0;
          for(size_t e=graph.in_row_ptr[v];e<graph.in_row_ptr[v+1];e++) {
            int u=graph.in_col_idx[e];
            incoming+=rank[u]/graph.out_degree[u];
          }
          next[v]=0.15/graph.n+0.85*incoming;
        }
      } else {
        #pragma omp parallel for schedule(guided,chunks[config])
        for(int v=0;v<graph.n;v++) {
          double incoming=0.0;
          for(size_t e=graph.in_row_ptr[v];e<graph.in_row_ptr[v+1];e++) {
            int u=graph.in_col_idx[e];
            incoming+=rank[u]/graph.out_degree[u];
          }
          next[v]=0.15/graph.n+0.85*incoming;
        }
      }
      double*tmp=rank;
      rank=next;
      next=tmp;
    }
    double elapsed=omp_get_wtime()-start;
    printf("%-8s %4d %.4f %.2f\n",names[config],chunks[config],elapsed,(graph.edge_count*15)/(elapsed*1e6));
  }
  free(rank);
  free(next);
  free(thread_time);
  free(thread_edges);
  free_graph(&graph);
  return 0;
}
