#!/bin/bash

#SBATCH -p Instruction
#SBATCH -J job11
#SBATCH -e job11.err
#SBATCH -o job11.out
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 00:05:00

export OMP_NUM_THREADS=4
./hello_omp_hb
