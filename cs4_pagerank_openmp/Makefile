CC=gcc
CFLAGS=-O3 -std=c11 -Wall -Wextra -fopenmp

all:serial parallel_pull parallel_push schedule_bench

serial:serial.c common.h
	$(CC) $(CFLAGS) -o serial serial.c

parallel_pull:parallel_pull.c common.h
	$(CC) $(CFLAGS) -o parallel_pull parallel_pull.c

parallel_push:parallel_push.c common.h
	$(CC) $(CFLAGS) -o parallel_push parallel_push.c

schedule_bench:schedule_bench.c common.h
	$(CC) $(CFLAGS) -o schedule_bench schedule_bench.c

clean:
	rm -f serial parallel_pull parallel_push schedule_bench
