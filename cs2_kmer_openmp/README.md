#Case Study 2:Parallel K-mer Counting

#Build
make

#Run
./kmer N k mode threads seed skewed

Modes:
serial
race
atomic
private
padded

Examples:
./kmer 100000000 8 serial 1 12345 0
./kmer 100000000 8 race 4 12345 0
./kmer 100000000 8 atomic 4 12345 0
./kmer 100000000 12 private 4 12345 0
./kmer 1000000 8 padded 4 12345 0

#Benchmarks
python3 benchmark.py 10000000 3 kmer_results.csv

#False-sharing check
./kmer 10000000 1 private 4 12345 0
./kmer 10000000 1 padded 4 12345 0
