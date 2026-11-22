#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdint.h>
#include<omp.h>
#include<time.h>

typedef enum {
  MODE_SERIAL,
  MODE_RACE,
  MODE_ATOMIC,
  MODE_PRIVATE,
  MODE_PADDED
}Mode;

static uint64_t now_ns(void) {
  struct timespec value;
  clock_gettime(CLOCK_MONOTONIC,&value);
  return (uint64_t)value.tv_sec*1000000000ULL+(uint64_t)value.tv_nsec;
}

static void encode_base(char*sequence,size_t length,uint64_t seed,int skewed) {
  uint64_t state=seed;
  for(size_t i=0;i<length;i++) {
    state=state*6364136223846793005ULL+1442695040888963407ULL;
    uint64_t random=state>>32;
    char base;
    if(skewed) {
      uint64_t choice=random%100;
      if(choice<60) {
        base='A';
      } else if(choice<80) {
        base='C';
      } else if(choice<90) {
        base='G';
      } else {
        base='T';
      }
    } else {
      base="ACGT"[random%4];
    }
    sequence[i]=base;
  }
}

static uint64_t encode_kmer(char*sequence,size_t start,int k) {
  uint64_t index=0;
  for(int i=0;i<k;i++) {
    char base=sequence[start+i];
    int code=(base=='C')?1:(base=='G')?2:(base=='T')?3:0;
    index=(index<<2)|code;
  }
  return index;
}

static uint64_t encode_first_kmer(char*sequence,int k) {
  return encode_kmer(sequence,0,k);
}

static int base_code(char base) {
  if(base=='C') {
    return 1;
  }
  if(base=='G') {
    return 2;
  }
  if(base=='T') {
    return 3;
  }
  return 0;
}

static uint64_t next_kmer(uint64_t index,char incoming,int k) {
  uint64_t mask=((uint64_t)1<<(2*k))-1;
  return ((index<<2)&mask)|base_code(incoming);
}

static Mode parse_mode(char*text) {
  if(!strcmp(text,"race")) {
    return MODE_RACE;
  }
  if(!strcmp(text,"atomic")) {
    return MODE_ATOMIC;
  }
  if(!strcmp(text,"private")) {
    return MODE_PRIVATE;
  }
  if(!strcmp(text,"padded")) {
    return MODE_PADDED;
  }
  return MODE_SERIAL;
}

static void count_serial(char*sequence,size_t length,int k,uint64_t*bins) {
  size_t windows=length-(size_t)k+1;
  uint64_t index=encode_first_kmer(sequence,k);
  bins[index]++;
  for(size_t i=1;i<windows;i++) {
    index=next_kmer(index,sequence[i+k-1],k);
    bins[index]++;
  }
}

static void count_race(char*sequence,size_t length,int k,uint64_t*bins) {
  size_t windows=length-(size_t)k+1;
  int threads=omp_get_max_threads();
  #pragma omp parallel
  {
    int thread=omp_get_thread_num();
    size_t begin=thread*windows/(size_t)threads;
    size_t end=(thread+1)*windows/(size_t)threads;
    if(begin<end) {
      uint64_t index=encode_kmer(sequence,begin,k);
      bins[index]++;
      for(size_t i=begin+1;i<end;i++) {
        index=next_kmer(index,sequence[i+k-1],k);
        bins[index]++;
      }
    }
  }
}

static void count_atomic(char*sequence,size_t length,int k,uint64_t*bins) {
  size_t windows=length-(size_t)k+1;
  int threads=omp_get_max_threads();
  #pragma omp parallel
  {
    int thread=omp_get_thread_num();
    size_t begin=thread*windows/(size_t)threads;
    size_t end=(thread+1)*windows/(size_t)threads;
    if(begin<end) {
      uint64_t index=encode_kmer(sequence,begin,k);
      #pragma omp atomic
      bins[index]++;
      for(size_t i=begin+1;i<end;i++) {
        index=next_kmer(index,sequence[i+k-1],k);
        #pragma omp atomic
        bins[index]++;
      }
    }
  }
}

static void count_private(char*sequence,size_t length,int k,uint64_t*bins,size_t bin_count,int padded) {
  int threads=omp_get_max_threads();
  size_t stride=bin_count;
  if(padded) {
    stride=((bin_count*sizeof(uint64_t)+63)/64)*64/sizeof(uint64_t);
  }
  uint64_t*local=(uint64_t*)calloc((size_t)threads*stride,sizeof(uint64_t));
  if(!local) {
    fprintf(stderr,"allocation failed\n");
    exit(1);
  }
  size_t windows=length-(size_t)k+1;
  #pragma omp parallel
  {
    int thread=omp_get_thread_num();
    uint64_t*hist=local+(size_t)thread*stride;
    size_t begin=thread*windows/(size_t)threads;
    size_t end=(thread+1)*windows/(size_t)threads;
    if(begin<end) {
      uint64_t index=encode_kmer(sequence,begin,k);
      hist[index]++;
      for(size_t i=begin+1;i<end;i++) {
        index=next_kmer(index,sequence[i+k-1],k);
        hist[index]++;
      }
    }
  }
  for(int thread=0;thread<threads;thread++) {
    uint64_t*hist=local+(size_t)thread*stride;
    for(size_t i=0;i<bin_count;i++) {
      bins[i]+=hist[i];
    }
  }
  free(local);
}

static int histogram_ok(uint64_t*bins,size_t bin_count,uint64_t expected) {
  uint64_t total=0;
  for(size_t i=0;i<bin_count;i++) {
    total+=bins[i];
  }
  return total==expected;
}

static void print_usage(void) {
  fprintf(stderr,"usage:kmer N k mode threads seed skewed\n");
  fprintf(stderr,"mode:serial race atomic private padded\n");
}

int main(int argc,char**argv) {
  size_t length=(argc>1)?(size_t)strtoull(argv[1],NULL,10):100000000ULL;
  int k=(argc>2)?atoi(argv[2]):8;
  Mode mode=(argc>3)?parse_mode(argv[3]):MODE_SERIAL;
  int threads=(argc>4)?atoi(argv[4]):4;
  uint64_t seed=(argc>5)?strtoull(argv[5],NULL,10):12345ULL;
  int skewed=(argc>6)?atoi(argv[6]):0;
  if(k<1||k>16||length<(size_t)k) {
    print_usage();
    return 1;
  }
  omp_set_num_threads(threads);
  char*sequence=(char*)malloc(length);
  if(!sequence) {
    fprintf(stderr,"allocation failed\n");
    return 1;
  }
  encode_base(sequence,length,seed,skewed);
  size_t bin_count=(size_t)1<<(2*k);
  uint64_t*bins=(uint64_t*)calloc(bin_count,sizeof(uint64_t));
  if(!bins) {
    fprintf(stderr,"allocation failed\n");
    free(sequence);
    return 1;
  }
  uint64_t start=now_ns();
  if(mode==MODE_SERIAL) {
    count_serial(sequence,length,k,bins);
  } else if(mode==MODE_RACE) {
    count_race(sequence,length,k,bins);
  } else if(mode==MODE_ATOMIC) {
    count_atomic(sequence,length,k,bins);
  } else if(mode==MODE_PRIVATE) {
    count_private(sequence,length,k,bins,bin_count,0);
  } else {
    count_private(sequence,length,k,bins,bin_count,1);
  }
  uint64_t finish=now_ns();
  uint64_t expected=length-(size_t)k+1;
  uint64_t total=0;
  for(size_t i=0;i<bin_count;i++) {
    total+=bins[i];
  }
  double seconds=(finish-start)/1000000000.0;
  printf("N=%zu k=%d mode=%s threads=%d skewed=%d time=%.9f total=%llu expected=%llu correct=%d throughput=%.3fM bins=%zu\n",
    length,k,argv[3]?argv[3]: "serial",threads,skewed,seconds,
    (unsigned long long)total,(unsigned long long)expected,histogram_ok(bins,bin_count,expected),
    (double)(length-(size_t)k+1)/1000000.0/seconds,bin_count);
  free(bins);
  free(sequence);
  return 0;
}
