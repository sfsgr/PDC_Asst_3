#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

int main(int argc,char**argv) {
  int scale=(argc>1)?atoi(argv[1]):16;
  size_t target=(argc>2)?strtoull(argv[2],NULL,10):1500000;
  int threads=(argc>3)?atoi(argv[3]):4;
  int schedule=(argc>4)?atoi(argv[4]):0;
  int chunk=(argc>5)?atoi(argv[5]):256;
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
  double damping=0.85;
  double tolerance=1e-7;
  int iterations=0;
  double residual=1.0;
  double start=omp_get_wtime();
  for(int iteration=1;iteration<=100;iteration++) {
    double dangling_sum=0.0;
    #pragma omp parallel for reduction(+:dangling_sum) schedule(static)
    for(int i=0;i<graph.dangling_count;i++) {
      dangling_sum+=rank[graph.dangling[i]];
    }
    double base=(1.0-damping)/graph.n+damping*dangling_sum/graph.n;
    residual=0.0;
    if(schedule==0) {
      #pragma omp parallel for reduction(+:residual) schedule(static,chunk)
      for(int v=0;v<graph.n;v++) {
        double incoming=0.0;
        for(size_t e=graph.in_row_ptr[v];e<graph.in_row_ptr[v+1];e++) {
          int u=graph.in_col_idx[e];
          incoming+=rank[u]/graph.out_degree[u];
        }
        next[v]=base+damping*incoming;
        residual+=fabs(next[v]-rank[v]);
      }
    } else if(schedule==1) {
      #pragma omp parallel for reduction(+:residual) schedule(dynamic,chunk)
      for(int v=0;v<graph.n;v++) {
        double incoming=0.0;
        for(size_t e=graph.in_row_ptr[v];e<graph.in_row_ptr[v+1];e++) {
          int u=graph.in_col_idx[e];
          incoming+=rank[u]/graph.out_degree[u];
        }
        next[v]=base+damping*incoming;
        residual+=fabs(next[v]-rank[v]);
      }
    } else {
      #pragma omp parallel for reduction(+:residual) schedule(guided,chunk)
      for(int v=0;v<graph.n;v++) {
        double incoming=0.0;
        for(size_t e=graph.in_row_ptr[v];e<graph.in_row_ptr[v+1];e++) {
          int u=graph.in_col_idx[e];
          incoming+=rank[u]/graph.out_degree[u];
        }
        next[v]=base+damping*incoming;
        residual+=fabs(next[v]-rank[v]);
      }
    }
    double*tmp=rank;
    rank=next;
    next=tmp;
    iterations=iteration;
    if(residual<tolerance) {
      break;
    }
  }
  double elapsed=omp_get_wtime()-start;
  printf("======================================================\n");
  printf("Case Study 4: Parallel Pull PageRank\n");
  printf("Vertices: %d Edges: %zu Threads: %d\n",graph.n,graph.edge_count,threads);
  printf("Schedule: %s Chunk: %d\n",schedule==0?"static":schedule==1?"dynamic":"guided",chunk);
  printf("Iterations: %d Final L1 Residual: %.8e\n",iterations,residual);
  printf("Compute Time: %.4f s Throughput: %.2f Million edges/s\n",elapsed,(graph.edge_count*iterations)/(elapsed*1e6));
  free(rank);
  free(next);
  free_graph(&graph);
  return 0;
}
