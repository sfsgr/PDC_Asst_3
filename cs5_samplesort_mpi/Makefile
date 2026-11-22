CC=gcc
MPICC=mpicc
CFLAGS=-O3 -std=c11 -Wall -Wextra

all:serial_sort parallel_sample_sort pingpong_bench

serial_sort:serial_sort.c common.h
	$(CC) $(CFLAGS) -o serial_sort serial_sort.c

parallel_sample_sort:parallel_sample_sort.c common.h
	$(MPICC) $(CFLAGS) -o parallel_sample_sort parallel_sample_sort.c

pingpong_bench:pingpong_bench.c
	$(MPICC) $(CFLAGS) -o pingpong_bench pingpong_bench.c

clean:
	rm -f serial_sort parallel_sample_sort pingpong_bench
