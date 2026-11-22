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

#Benchmark
python3 benchmark.py 100000000 3 sample_sort_results.csv

#Measured 100-million-key results
See sample_sort_100m_results.csv.

At p=4, oversample=1:
Uniform: max/mean=1.001
Skewed: max/mean=1.000

Uniform at p=4, oversample=4:
Min=23076323
Max=30760328
Max/mean=1.230

Measured ping-pong parameters:
Alpha=0.388us
Beta=0.1931ns/byte
Bandwidth=4.82GB/s
