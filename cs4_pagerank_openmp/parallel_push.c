#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

int main(int argc,char**argv) {
  int scale=(argc>1)?atoi(argv[1]):16;
  size_t target=(argc>2)?strtoull(argv[2],NULL,10):1500000;
  int threads=(argc>3)?atoi(argv[3]):4;
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
    #pragma omp parallel for schedule(static)
    for(int v=0;v<graph.n;v++) {
      next[v]=base;
    }
    #pragma omp parallel for schedule(static)
    for(int u=0;u<graph.n;u++) {
      if(graph.out_degree[u]==0) {
        continue;
      }
      double contribution=damping*rank[u]/graph.out_degree[u];
      for(size_t e=graph.row_ptr[u];e<graph.row_ptr[u+1];e++) {
        int v=graph.col_idx[e];
        #pragma omp atomic
        next[v]+=contribution;
      }
    }
    residual=0.0;
    #pragma omp parallel for reduction(+:residual) schedule(static)
    for(int v=0;v<graph.n;v++) {
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
  double elapsed=omp_get_wtime()-start;
  printf("======================================================\n");
  printf("Case Study 4: Parallel Push PageRank\n");
  printf("Vertices: %d Edges: %zu Threads: %d\n",graph.n,graph.edge_count,threads);
  printf("Iterations: %d Final L1 Residual: %.8e\n",iterations,residual);
  printf("Compute Time: %.4f s Throughput: %.2f Million edges/s\n",elapsed,(graph.edge_count*iterations)/(elapsed*1e6));
  double rank_sum=0.0;
  for(int i=0;i<graph.n;i++) {
    rank_sum+=rank[i];
  }
  printf("Rank Sum: %.10f (tolerance=1e-7)\n",rank_sum);
  free(rank);
  free(next);
  free_graph(&graph);
  return 0;
}
