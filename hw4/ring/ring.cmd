#PBS -S /bin/bash
#PBS -q newest
#PBS -N ring
#PBS -l nodes=1:ppn=5
#PBS -l walltime=00:05:00

cd $PBS_O_WORKDIR
mpirun -np 5 ring
