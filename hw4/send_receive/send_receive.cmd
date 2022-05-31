#PBS -S /bin/bash
#PBS -q newest
#PBS -N send_receive
#PBS -l nodes=1:ppn=2
#PBS -l walltime=00:05:00

cd $PBS_O_WORKDIR
mpirun -np 2 send_receive
