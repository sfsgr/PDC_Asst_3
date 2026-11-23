
#include "mc.h"
#include <stdlib.h>
#include <time.h>

int main(int argc,char**argv){
    unsigned long long samples=100000000;
    unsigned long long seed=12345;
    if(argc>1)samples=strtoull(argv[1],NULL,10);
    if(argc>2)seed=strtoull(argv[2],NULL,10);
    if(samples<2){
        printf("Samples must be at least 2\n");
        return 1;
    }
    RNG rng;
    rng_init(&rng,seed,0);
    double sum=0.0;
    double sum2=0.0;
    clock_t start=clock();
    for(unsigned long long i=0;i<samples;i++){
        double payoff=get_payoff(&rng);
        sum+=payoff;
        sum2+=payoff*payoff;
    }
    double elapsed=(double)(clock()-start)/CLOCKS_PER_SEC;
    print_result("serial",samples,seed,1,sum,sum2,elapsed);
    return 0;
}
