#define _POSIX_C_SOURCE 200809L
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

typedef struct {
  int n;
  int steps;
  int k;
  double alpha;
}Params;

static double*allocate_grid(int n) {
  double*p=(double*)malloc((size_t)n*n*sizeof(double));
  if(!p) {
    fprintf(stderr,"allocation failed\n");
    exit(1);
  }
  return p;
}

static void initialize_grid(double*a,int n) {
  memset(a,0,(size_t)n*n*sizeof(double));
  for(int j=0;j<n;j++) {
    a[j]=100.0;
    a[(n-1)*n+j]=100.0;
  }
  for(int i=0;i<n;i++) {
    a[i*n]=100.0;
    a[i*n+n-1]=100.0;
  }
}

static void copy_boundaries(double*a,double*b,int n) {
  for(int j=0;j<n;j++) {
    b[j]=a[j];
    b[(n-1)*n+j]=a[(n-1)*n+j];
  }
  for(int i=0;i<n;i++) {
    b[i*n]=a[i*n];
    b[i*n+n-1]=a[i*n+n-1];
  }
}

static double run_step(double*a,double*b,Params*params) {
  int n=params->n;
  double alpha=params->alpha;
  double max_change=0.0;
  for(int i=1;i<n-1;i++) {
    for(int j=1;j<n-1;j++) {
      double old=a[i*n+j];
      double next=old+alpha*(a[(i-1)*n+j]+a[(i+1)*n+j]+a[i*n+j-1]+a[i*n+j+1]-4.0*old);
      b[i*n+j]=next;
      double change=fabs(next-old);
      if(change>max_change) {
        max_change=change;
      }
    }
  }
  copy_boundaries(a,b,n);
  return max_change;
}

static int parse_params(int argc,char**argv,Params*params) {
  params->n=(argc>1)?atoi(argv[1]):1024;
  params->steps=(argc>2)?atoi(argv[2]):500;
  params->alpha=(argc>3)?atof(argv[3]):0.1;
  params->k=(argc>4)?atoi(argv[4]):50;
  if(params->n<3||params->steps<0||params->k<1||params->alpha<0.0||params->alpha>0.25) {
    fprintf(stderr,"usage:serial N T alpha K\n");
    return 0;
  }
  return 1;
}

int main(int argc,char**argv) {
  Params params;
  if(!parse_params(argc,argv,&params)) {
    return 1;
  }
  double*a=allocate_grid(params.n);
  double*b=allocate_grid(params.n);
  initialize_grid(a,params.n);
  memcpy(b,a,(size_t)params.n*params.n*sizeof(double));
  double last_metric=0.0;
  struct timespec start,end;
  clock_gettime(CLOCK_MONOTONIC,&start);
  for(int step=0;step<params.steps;step++) {
    double change=run_step(a,b,&params);
    double*tmp=a;
    a=b;
    b=tmp;
    if((step+1)%params.k==0) {
      last_metric=change;
    }
  }
  clock_gettime(CLOCK_MONOTONIC,&end);
  double seconds=(double)(end.tv_sec-start.tv_sec)+(double)(end.tv_nsec-start.tv_nsec)*1e-9;
  double sum=0.0;
  for(int i=0;i<params.n*params.n;i++) {
    sum+=a[i];
  }
  printf("N=%d T=%d alpha=%.17g K=%d time=%.9f sum=%.17g max_change=%.17g\n",
    params.n,params.steps,params.alpha,params.k,seconds,sum,last_metric);
  free(a);
  free(b);
  return 0;
}
