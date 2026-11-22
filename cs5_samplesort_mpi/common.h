#ifndef COMMON_H
#define COMMON_H

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>

static uint64_t sample_rng(uint64_t*state) {
  uint64_t z=(*state+=0x9e3779b97f4a7c15ULL);
  z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;
  z=(z^(z>>27))*0x94d049bb133111ebULL;
  return z^(z>>31);
}

static void generate_uniform_data(uint32_t*data,size_t count,uint64_t seed) {
  uint64_t state=seed;
  for(size_t i=0;i<count;i++) {
    data[i]=(uint32_t)sample_rng(&state);
  }
}

static void generate_skewed_data(uint32_t*data,size_t count,uint64_t seed) {
  uint64_t state=seed;
  for(size_t i=0;i<count;i++) {
    uint32_t value=(uint32_t)sample_rng(&state);
    if((sample_rng(&state)%100)<80) {
      data[i]=(value%(1U<<24));
    } else {
      data[i]=value;
    }
  }
}

static int uint32_compare(const void*a,const void*b) {
  uint32_t x=*(const uint32_t*)a;
  uint32_t y=*(const uint32_t*)b;
  return (x>y)-(x<y);
}

static int verify_local_sorted(const uint32_t*data,size_t count) {
  for(size_t i=1;i<count;i++) {
    if(data[i-1]>data[i]) {
      return 0;
    }
  }
  return 1;
}

#endif
