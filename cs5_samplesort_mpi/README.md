#Case Study 5:Sample Sort

#Build
make

#Serial
./serial_sort N distribution
./serial_sort 10000000 0

#MPI
mpirun -n p ./parallel_sample_sort N distribution oversample
mpirun -n 4 ./parallel_sample_sort 10000000 0 1
mpirun -n 4 ./parallel_sample_sort 10000000 1 1

#Ping-pong
mpirun -n 2 ./pingpong_bench
