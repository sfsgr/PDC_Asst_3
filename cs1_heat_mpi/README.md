#Case Study 1:Heat Diffusion with MPI

#Build
sudo apt update
sudo apt install build-essential openmpi-bin libopenmpi-dev
make

#Serial
./serial N T alpha K
Example:
./serial 1024 500 0.1 50

#MPI
mpirun -n p ./parallel N T alpha K 1d|2d blocking|nonblocking
Examples:
mpirun -n 4 ./parallel 1024 500 0.1 50 1d blocking
mpirun -n 4 ./parallel 1024 500 0.1 50 2d nonblocking

The fixed boundary value is 100.0 and the interior starts at 0.0. Parallel floating-point results should match the serial sum within 1e-9. K controls the MPI_Allreduce convergence check. K=1 checks every step and K=50 checks every 50 steps.

#Plot scaling
Install the plot dependency:
sudo apt install python3-matplotlib

python3 plot_scaling.py strong results.csv strong_scaling.png
python3 plot_scaling.py weak weak_results.csv weak_scaling.png

CSV columns for strong scaling: ranks,t1,time,n
CSV columns for weak scaling: ranks,time,n
