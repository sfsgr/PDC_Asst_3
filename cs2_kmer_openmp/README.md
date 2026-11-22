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

#Critical analysis

##Correctness
The race version at k=8 and 4 threads counted 99,946,568 k-mers instead of 99,999,993.
The missing 53,425 increments are lost updates. Two threads can both read the same histogram counter, then both write the same incremented value.

Chunk boundaries are handled by assigning whole k-mer windows to threads and starting each chunk with encode_kmer at its first window.
The remaining windows use a rolling 2-bit index, so no k-mer is duplicated or omitted at a chunk boundary.

##Atomic versus private
For 100 million characters and 4 threads:

|k|Atomic|Private|Padded private|
|---|---|---|---|
|8|61.966M/s|237.719M/s|248.884M/s|
|12|51.713M/s|24.718M/s|Not measured|

Private is faster for k=8 because 65,536 counters are small enough to remain cache-resident per thread.
Atomic is faster for k=12 because private storage is 16,777,216 counters per thread and exceeds cache capacity.
The private histogram memory is 4^k*threads*8 bytes.

|k|One private histogram|Extra at 4 threads|Extra at 8 threads|
|---|---|---|---|
|8|512KiB|2MiB|4MiB|
|12|128MiB|512MiB|1GiB|

##False sharing
At k=1 and 2 threads, private histograms used 0.107876888 seconds and padded histograms used 0.071587092 seconds.
Padding improved time by about 50.7 percent and throughput from 92.698M/s to 139.690M/s.
This is false sharing: adjacent per-thread rows share cache lines and create write traffic between cores.

##Memory or compute
The rolling implementation reads one sequence byte for each k-mer and reads/writes one 8-byte histogram counter.
A simple 100-million-character estimate is 9 bytes read per k-mer plus histogram write traffic.
Measured k=8 throughput is 248.884M k-mers/s for padded private at 4 threads.
Scaling flattens at 8 threads because of cache pressure, merge cost, and oversubscription beyond 4 physical cores.
