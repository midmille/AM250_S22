#PBS -S /bin/bash
#PBS -q newest 
#PBS -N hello_world
#PBS -l nodes=1:ppn=6
#PBS -l walltime=00:05:00

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=6
./hello_world
