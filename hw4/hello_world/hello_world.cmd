#PBS -S /bin/bash
#PBS -q newest 
#PBS -N hello_world
#PBS -l nodes=1:ppn=4
#PBS -l walltime=00:05:00

cd $PBS_O_WORKDIR
mpirun -np 4 hello_world
