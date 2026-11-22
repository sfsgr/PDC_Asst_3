#Case Study 4:PageRank

#Build
make

#Run
./serial N scale target_edges
./serial 16 1500000
./parallel_pull 16 1500000 4 0 256
./parallel_push 16 1500000 4
./schedule_bench 16 1500000 8
