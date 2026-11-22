CC=gcc
CFLAGS=-O3 -std=c11 -Wall -Wextra -fopenmp

.PHONY:all clean

all:kmer

kmer:kmer.c
	$(CC) $(CFLAGS) -o kmer kmer.c

clean:
	rm -f kmer
