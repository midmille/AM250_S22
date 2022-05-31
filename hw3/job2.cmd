#PBS -S /bin/bash
#PBS -q newest
#PBS -N job2
#PBS -l nodes=2:ppn=4
#PBS -l walltime=00:05:00

cd $PBS_O_WORKDIR 
mpirun -np 8 hello_mpi_grape
