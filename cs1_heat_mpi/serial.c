#define _POSIX_C_SOURCE 200809L
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

static double*allocate_grid(int n){
double*p=(double*)malloc((size_t)n*n*sizeof(double));
if(!p){fprintf(stderr,"allocation failed\n");exit(1);}
return p;
}

static void initialize(double*a,int n){
memset(a,0,(size_t)n*n*sizeof(double));
for(int j=0;j<n;j++){
a[j]=100.0;
a[(n-1)*n+j]=100.0;
}
for(int i=0;i<n;i++){
a[i*n]=100.0;
a[i*n+n-1]=100.0;
}
}

int main(int argc,char**argv){
int n=1024;
int steps=500;
int k=50;
double alpha=0.1;
if(argc>1)n=atoi(argv[1]);
if(argc>2)steps=atoi(argv[2]);
if(argc>3)alpha=atof(argv[3]);
if(argc>4)k=atoi(argv[4]);
if(n<3||steps<0||k<1||alpha<0.0||alpha>0.25){fprintf(stderr,"usage: serial N T alpha K\n");return 1;}
double*a=allocate_grid(n);
double*b=allocate_grid(n);
initialize(a,n);
memcpy(b,a,(size_t)n*n*sizeof(double));
double last_metric=0.0;
struct timespec t0,t1;
clock_gettime(CLOCK_MONOTONIC,&t0);
for(int step=0;step<steps;step++){
double max_change=0.0;
for(int i=1;i<n-1;i++){
for(int j=1;j<n-1;j++){
double old=a[i*n+j];
double next=old+alpha*(a[(i-1)*n+j]+a[(i+1)*n+j]+a[i*n+j-1]+a[i*n+j+1]-4.0*old);
b[i*n+j]=next;
double change=fabs(next-old);
if(change>max_change)max_change=change;
}
}
for(int j=0;j<n;j++){
b[j]=a[j];
b[(n-1)*n+j]=a[(n-1)*n+j];
}
for(int i=0;i<n;i++){
b[i*n]=a[i*n];
b[i*n+n-1]=a[i*n+n-1];
}
double*tmp=a;a=b;b=tmp;
if((step+1)%k==0)last_metric=max_change;
}
clock_gettime(CLOCK_MONOTONIC,&t1);
double seconds=(double)(t1.tv_sec-t0.tv_sec)+(double)(t1.tv_nsec-t0.tv_nsec)*1e-9;
double sum=0.0;
for(int i=0;i<n*n;i++)sum+=a[i];
printf("N=%d T=%d alpha=%.17g K=%d time=%.9f sum=%.17g max_change=%.17g\n",n,steps,alpha,k,seconds,sum,last_metric);
free(a);
free(b);
return 0;
}
