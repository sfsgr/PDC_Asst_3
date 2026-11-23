
#include "mc.h"
#include <stdlib.h>
#include <string.h>
#include <omp.h>

int main(int argc,char**argv){
    unsigned long long samples=100000000;
    unsigned long long seed=12345;
    int threads=omp_get_max_threads();
    int timed=0;
    int bad=0;
    double seconds=10.0;

    if(argc>1&&strcmp(argv[1],"time")==0){
        if(argc!=5){
            printf("Usage: ./parallel time seconds seed threads\n");
            return 1;
        }
        timed=1;
        seconds=atof(argv[2]);
        seed=strtoull(argv[3],NULL,10);
        threads=atoi(argv[4]);
    }else if(argc>1&&strcmp(argv[1],"same-seed")==0){
        if(argc!=5){
            printf("Usage: ./parallel same-seed M seed threads\n");
            return 1;
        }
        bad=1;
        samples=strtoull(argv[2],NULL,10);
        seed=strtoull(argv[3],NULL,10);
        threads=atoi(argv[4]);
    }else{
        if(argc>1)samples=strtoull(argv[1],NULL,10);
        if(argc>2)seed=strtoull(argv[2],NULL,10);
        if(argc>3)threads=atoi(argv[3]);
    }

    if(threads<1||(!timed&&samples<2)||(timed&&seconds<=0.0)){
        printf("Invalid parameters\n");
        return 1;
    }

    omp_set_dynamic(0);
    omp_set_num_threads(threads);

    double*sums=calloc(threads,sizeof(double));
    double*squares=calloc(threads,sizeof(double));
    unsigned long long*counts=calloc(threads,sizeof(unsigned long long));

    if(!sums||!squares||!counts){
        printf("Memory allocation failed\n");
        free(sums);
        free(squares);
        free(counts);
        return 1;
    }

    int actual_threads=0;
    double start=omp_get_wtime();

    #pragma omp parallel
    {
        int tid=omp_get_thread_num();
        int workers=omp_get_num_threads();
        RNG rng;
        rng_init(&rng,seed,bad?0:(uint64_t)tid);

        double local_sum=0.0;
        double local_sum2=0.0;
        unsigned long long local_count=0;

        if(timed){
            while(1){
                if((local_count&4095ULL)==0&&omp_get_wtime()-start>=seconds)break;
                double payoff=get_payoff(&rng);
                local_sum+=payoff;
                local_sum2+=payoff*payoff;
                local_count++;
            }
        }else{
            unsigned long long target=samples/workers+(tid<samples%workers);
            for(unsigned long long i=0;i<target;i++){
                double payoff=get_payoff(&rng);
                local_sum+=payoff;
                local_sum2+=payoff*payoff;
                local_count++;
            }
        }

        sums[tid]=local_sum;
        squares[tid]=local_sum2;
        counts[tid]=local_count;

        #pragma omp single
        actual_threads=workers;
    }

    double sum=0.0;
    double sum2=0.0;
    unsigned long long total=0;

    for(int i=0;i<actual_threads;i++){
        sum+=sums[i];
        sum2+=squares[i];
        total+=counts[i];
    }

    double elapsed=omp_get_wtime()-start;
    const char*mode=timed?"timed":(bad?"same-seed":"parallel");

    print_result(mode,total,seed,actual_threads,sum,sum2,elapsed);

    free(sums);
    free(squares);
    free(counts);
    return 0;
}
