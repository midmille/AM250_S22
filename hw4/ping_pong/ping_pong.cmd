#PBS -S /bin/bash
#PBS -q newest
#PBS -N ping_pong
#PBS -l nodes=1:ppn=2
#PBS -l walltime=00:05:00

cd $PBS_O_WORKDIR
mpirun -np 2 ping_pong
