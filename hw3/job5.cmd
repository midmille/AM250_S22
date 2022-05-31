#!/bin/bash

#SBATCH -p Instruction
#SBATCH -J job5
#SBATCH -e job5.err
#SBATCH -o job5.out
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 00:05:00

mpirun -np 4 hello_mpi_hb
