CC=mpicc
.CFLAGS=-O3 -std=c11 -Wall -Wextra

.PHONY:all clean

all:serial parallel

serial:serial.c
 gcc $(CFLAGS) -o serial serial.c -lm

parallel:parallel.c
 $(CC) $(CFLAGS) -o parallel parallel.c -lm

clean:
 rm -f serial parallel
