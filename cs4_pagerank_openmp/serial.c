#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

static double elapsed(struct timespec start,struct timespec end) {
  return (double)(end.tv_sec-start.tv_sec)+(double)(end.tv_nsec-start.tv_nsec)/1e9;
}

int main(int argc,char**argv) {
  int scale=(argc>1)?atoi(argv[1]):16;
  size_t target=(argc>2)?strtoull(argv[2],NULL,10):1500000;
  double damping=0.85;
  double tolerance=1e-7;
  int max_iterations=100;
  struct timespec start,end;
  clock_gettime(CLOCK_MONOTONIC,&start);
  CSRGraph graph=generate_rmat_graph(scale,target,42ULL);
  clock_gettime(CLOCK_MONOTONIC,&end);
  double generation=elapsed(start,end);
  double*rank=(double*)malloc((size_t)graph.n*sizeof(double));
  double*next=(double*)malloc((size_t)graph.n*sizeof(double));
  if(!rank||!next) {
    fprintf(stderr,"allocation failed\n");
    return 1;
  }
  for(int i=0;i<graph.n;i++) {
    rank[i]=1.0/graph.n;
  }
  clock_gettime(CLOCK_MONOTONIC,&start);
  int iterations=0;
  double residual=1.0;
  for(int iteration=1;iteration<=max_iterations;iteration++) {
    double dangling_sum=0.0;
    for(int i=0;i<graph.dangling_count;i++) {
      dangling_sum+=rank[graph.dangling[i]];
    }
    double base=(1.0-damping)/graph.n+damping*dangling_sum/graph.n;
    residual=0.0;
    for(int v=0;v<graph.n;v++) {
      double incoming=0.0;
      for(size_t e=graph.in_row_ptr[v];e<graph.in_row_ptr[v+1];e++) {
        int u=graph.in_col_idx[e];
        incoming+=rank[u]/graph.out_degree[u];
      }
      next[v]=base+damping*incoming;
      residual+=fabs(next[v]-rank[v]);
    }
    double*tmp=rank;
    rank=next;
    next=tmp;
    iterations=iteration;
    if(residual<tolerance) {
      break;
    }
  }
  clock_gettime(CLOCK_MONOTONIC,&end);
  double compute=elapsed(start,end);
  printf("======================================================\n");
  printf("Case Study 4: Serial PageRank\n");
  printf("Vertices: %d\n",graph.n);
  printf("Edges: %zu\n",graph.edge_count);
  printf("Generation Time: %.4f s\n",generation);
  printf("Iterations: %d\n",iterations);
  printf("Final L1 Residual: %.8e\n",residual);
  printf("Compute Time: %.4f s\n",compute);
  printf("Throughput: %.2f Million edges/s\n",(graph.edge_count*iterations)/(compute*1e6));
  free(rank);
  free(next);
  free_graph(&graph);
  return 0;
}
