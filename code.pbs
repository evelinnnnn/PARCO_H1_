#!/bin/bash
# Job name
#PBS -N test
#PBS -j oe
# Queue name
#PBS -q short_cpuQ
# Set the maximum wall time
#PBS -l walltime=0:10:00
# Number of nodes, cpus, mpi processors and amount of memory
#PBS -l select=1:ncpus=64:ompthreads=64:mem=1mb

# Modules for C
module load gcc91
gcc() {
    gcc-9.1.0 "$@"
}
gcc --version
# Select the working directory
cd /home/evelin.begher/H1

# The code should be compiled before submitting the job
gcc -o exeC H1.c -fopenmp

# Loop over different thread and matrix dimention counts
for n in 16 32 64 128 256 512 1024 2048 4096; do

    echo >> H1_output.log
    echo "                $n x $n           " >> H1_output.log
    echo  >> H1_output.log
 # Set the number of threads and run the program
for threads in 2 4 8 16; do
    echo "Running with $threads threads" >> H1_output.log
      for i in {1..10}; do
      echo "::::::::::::::::::::::::::::  $i  ::::::::::::::::::::::::::::::::::::::::"  >> H1_output.log
      OMP_NUM_THREADS=$threads ./exeC $n >> H1_output.log
      
      done 
        echo >> H1_output.log
        echo  >> H1_output.log
  done
done

