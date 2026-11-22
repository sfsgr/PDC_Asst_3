#ifndef COMMON_H
#define COMMON_H

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct {
  int n;
  size_t edge_count;
  size_t*row_ptr;
  size_t*in_row_ptr;
  int*col_idx;
  int*in_col_idx;
  int*out_degree;
  int*in_degree;
  int*dangling;
  int dangling_count;
}CSRGraph;

static uint64_t graph_rng(uint64_t*state) {
  uint64_t z=(*state+=0x9e3779b97f4a7c15ULL);
  z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;
  z=(z^(z>>27))*0x94d049bb133111ebULL;
  return z^(z>>31);
}

static double graph_uniform(uint64_t*state) {
  return (graph_rng(state)>>11)*(1.0/9007199254740992.0);
}

typedef struct {
  int src;
  int dst;
}Edge;

static int edge_compare(const void*a,const void*b) {
  const Edge*x=(const Edge*)a;
  const Edge*y=(const Edge*)b;
  if(x->src!=y->src) {
    return (x->src>y->src)-(x->src<y->src);
  }
  return (x->dst>y->dst)-(x->dst<y->dst);
}

static CSRGraph generate_rmat_graph(int scale,size_t target_edges,uint64_t seed) {
  CSRGraph g;
  g.n=1<<scale;
  g.edge_count=0;
  g.row_ptr=(size_t*)calloc((size_t)g.n+1,sizeof(size_t));
  g.in_row_ptr=(size_t*)calloc((size_t)g.n+1,sizeof(size_t));
  g.col_idx=NULL;
  g.in_col_idx=NULL;
  g.out_degree=(int*)calloc((size_t)g.n,sizeof(int));
  g.in_degree=(int*)calloc((size_t)g.n,sizeof(int));
  g.dangling=(int*)malloc((size_t)g.n*sizeof(int));
  g.dangling_count=0;
  Edge*edges=(Edge*)malloc(target_edges*sizeof(Edge));
  if(!g.row_ptr||!g.in_row_ptr||!g.out_degree||!g.in_degree||!g.dangling||!edges) {
    fprintf(stderr,"allocation failed\n");
    exit(1);
  }
  uint64_t state=seed;
  size_t edge_count=0;
  for(size_t e=0;e<target_edges;e++) {
    int u=0;
    int v=0;
    for(int depth=0;depth<scale;depth++) {
      double p=graph_uniform(&state);
      int bit=1<<(scale-depth-1);
      if(p<0.57) {
      } else if(p<0.76) {
        v+=bit;
      } else if(p<0.95) {
        u+=bit;
      } else {
        u+=bit;
        v+=bit;
      }
    }
    if(u!=v) {
      edges[edge_count].src=u;
      edges[edge_count].dst=v;
      edge_count++;
    }
  }
  qsort(edges,edge_count,sizeof(Edge),edge_compare);
  size_t unique_count=0;
  for(size_t i=0;i<edge_count;i++) {
    if(i==0||edges[i].src!=edges[unique_count-1].src||edges[i].dst!=edges[unique_count-1].dst) {
      edges[unique_count++]=edges[i];
    }
  }
  g.edge_count=unique_count;
  for(size_t i=0;i<unique_count;i++) {
    g.out_degree[edges[i].src]++;
  }
  for(int i=0;i<g.n;i++) {
    g.row_ptr[i+1]=g.row_ptr[i]+(size_t)g.out_degree[i];
  }
  g.col_idx=(int*)malloc(unique_count*sizeof(int));
  size_t*cursor=(size_t*)malloc((size_t)(g.n+1)*sizeof(size_t));
  for(int i=0;i<=g.n;i++) {
    cursor[i]=g.row_ptr[i];
  }
  for(size_t i=0;i<unique_count;i++) {
    g.col_idx[cursor[edges[i].src]++]=edges[i].dst;
  }
  free(cursor);
  for(size_t i=0;i<unique_count;i++) {
    g.in_degree[edges[i].dst]++;
  }
  for(int i=0;i<g.n;i++) {
    g.in_row_ptr[i+1]=g.in_row_ptr[i]+(size_t)g.in_degree[i];
  }
  g.in_col_idx=(int*)malloc(unique_count*sizeof(int));
  cursor=(size_t*)malloc((size_t)(g.n+1)*sizeof(size_t));
  for(int i=0;i<=g.n;i++) {
    cursor[i]=g.in_row_ptr[i];
  }
  for(size_t i=0;i<unique_count;i++) {
    g.in_col_idx[cursor[edges[i].dst]++]=edges[i].src;
  }
  free(cursor);
  free(edges);
  for(int i=0;i<g.n;i++) {
    if(g.out_degree[i]==0) {
      g.dangling[g.dangling_count++]=i;
    }
  }
  return g;
}

static void free_graph(CSRGraph*g) {
  free(g->row_ptr);
  free(g->in_row_ptr);
  free(g->col_idx);
  free(g->in_col_idx);
  free(g->out_degree);
  free(g->in_degree);
  free(g->dangling);
}

#endif
