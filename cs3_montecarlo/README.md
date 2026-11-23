
# Case Study 3 - Monte Carlo Option Pricing

## Model
OpenMP

## Parameters
S0=100
K=100
r=0.05
sigma=0.20
T=1
M=100000000
Seed=12345

## Compile
make clean
make

## Serial
./serial 100000000 12345

## Parallel
./parallel 100000000 12345 1
./parallel 100000000 12345 2
./parallel 100000000 12345 4
./parallel 100000000 12345 8

## Same Seed Experiment
./parallel same-seed 100000000 12345 8

## Fixed Time Experiment
./parallel time 10 12345 1
./parallel time 10 12345 8

## Clean
make clean
