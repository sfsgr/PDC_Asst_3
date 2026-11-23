
CC=gcc
CFLAGS=-O3 -std=c11 -Wall -Wextra
OMP=-fopenmp

all: serial parallel

serial: serial.c mc.h
	$(CC) $(CFLAGS) -o serial serial.c -lm

parallel: parallel.c mc.h
	$(CC) $(CFLAGS) $(OMP) -o parallel parallel.c -lm

clean:
	rm -f serial parallel
