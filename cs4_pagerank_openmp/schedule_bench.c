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
  free_graph(&graph);
  return 0;
}
