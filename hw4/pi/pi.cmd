#PBS -S /bin/bash
#PBS -q newest
#PBS -N pi
#PBS -l nodes=1:ppn=10
#PBS -l walltime=00:10:00

cd $PBS_O_WORKDIR
mpirun -np 10 pi
