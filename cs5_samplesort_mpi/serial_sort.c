#define _POSIX_C_SOURCE 200809L
#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double elapsed_seconds(struct timespec start,struct timespec end) {
  return (double)(end.tv_sec-start.tv_sec)+(double)(end.tv_nsec-start.tv_nsec)/1e9;
}

int main(int argc,char**argv) {
  size_t n=(argc>1)?strtoull(argv[1],NULL,10):10000000;
  int distribution=(argc>2)?atoi(argv[2]):0;
  uint32_t*data=(uint32_t*)malloc(n*sizeof(uint32_t));
  if(!data) {
    fprintf(stderr,"allocation failed\n");
    return 1;
  }
  struct timespec generation_start,generation_end;
  clock_gettime(CLOCK_MONOTONIC,&generation_start);
  if(distribution==0) {
    generate_uniform_data(data,n,123456789ULL);
  } else {
    generate_skewed_data(data,n,123456789ULL);
  }
  clock_gettime(CLOCK_MONOTONIC,&generation_end);
  double generation_time=elapsed_seconds(generation_start,generation_end);
  struct timespec sort_start,sort_end;
  clock_gettime(CLOCK_MONOTONIC,&sort_start);
  qsort(data,n,sizeof(uint32_t),uint32_compare);
  clock_gettime(CLOCK_MONOTONIC,&sort_end);
  double sort_time=elapsed_seconds(sort_start,sort_end);
  int sorted=verify_local_sorted(data,n);
  printf("======================================================\n");
  printf("Case Study 5: Serial Baseline (qsort)\n");
  printf("======================================================\n");
  printf("Total Keys (N): %zu\n",n);
  printf("Distribution: %s\n",distribution==0?"Uniform Random":"Skewed");
  printf("Generation Time: %.4f s\n",generation_time);
  printf("Sorting Time: %.4f s\n",sort_time);
  printf("Throughput: %.2f Million keys/s\n",(n/1e6)/sort_time);
  printf("Sortedness Check: %s\n",sorted?"PASSED":"FAILED");
  printf("======================================================\n");
  free(data);
  return sorted?0:1;
}
