#!/bin/bash

#SBATCH -p Instruction
#SBATCH -J job6
#SBATCH -e job6.err
#SBATCH -o job6.out
#SBATCH -N 2
#SBATCH -n 8
#SBATCH -t 00:05:00

mpirun -np 8 hello_mpi_hb
