#!/bin/bash

#SBATCH -p Instruction
#SBATCH -J job12
#SBATCH -e job12.err
#SBATCH -o job12.out
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 00:05:00

export OMP_NUM_THREADS=8
./hello_omp_hb
