
#ifndef MC_H
#define MC_H
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#define S0 100.0
#define K 100.0
#define R 0.05
#define SIGMA 0.20
#define T 1.0
#define PI 3.14159265358979323846

typedef struct{
    uint64_t state;
    uint64_t inc;
}RNG;

static uint64_t mix64(uint64_t x){
    x^=x>>30;
    x*=0xbf58476d1ce4e5b9ULL;
    x^=x>>27;
    x*=0x94d049bb133111ebULL;
    x^=x>>31;
    return x;
}

static uint32_t pcg32(RNG*rng){
    uint64_t old=rng->state;
    rng->state=old*6364136223846793005ULL+rng->inc;
    uint32_t x=(uint32_t)(((old>>18)^old)>>27);
    uint32_t rot=(uint32_t)(old>>59);
    return (x>>rot)|(x<<((-rot)&31));
}

static void rng_init(RNG*rng,uint64_t seed,uint64_t worker){
    uint64_t x=seed+0x9e3779b97f4a7c15ULL*(worker+1);
    rng->state=0;
    rng->inc=(mix64(x+0xda3e39cb94b95bdbULL)<<1)|1;
    pcg32(rng);
    rng->state+=mix64(x);
    pcg32(rng);
}

static double uniform(RNG*rng){
    return ((double)pcg32(rng)+0.5)/4294967296.0;
}

static double get_payoff(RNG*rng){
    double u1=uniform(rng);
    double u2=uniform(rng);
    double z=sqrt(-2.0*log(u1))*cos(2.0*PI*u2);
    double st=S0*exp((R-SIGMA*SIGMA/2.0)*T+SIGMA*sqrt(T)*z);
    return exp(-R*T)*fmax(st-K,0.0);
}

static double black_scholes(void){
    double d1=(log(S0/K)+(R+SIGMA*SIGMA/2.0)*T)/(SIGMA*sqrt(T));
    double d2=d1-SIGMA*sqrt(T);
    double n1=0.5*erfc(-d1/sqrt(2.0));
    double n2=0.5*erfc(-d2/sqrt(2.0));
    return S0*n1-K*exp(-R*T)*n2;
}

static void print_result(const char*mode,unsigned long long samples,unsigned long long seed,int threads,double sum,double sum2,double elapsed){
    double n=(double)samples;
    double price=sum/n;
    double variance=(sum2-sum*sum/n)/(n-1.0);
    if(variance<0.0)variance=0.0;
    double se=sqrt(variance/n);
    double bs=black_scholes();
    printf("mode=%s samples=%llu seed=%llu threads=%d price=%.10f se=%.10f bs=%.10f error=%.10f time=%.6f rate=%.2f\n",mode,samples,seed,threads,price,se,bs,fabs(price-bs),elapsed,samples/elapsed);
}
#endif
