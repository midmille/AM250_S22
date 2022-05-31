#PBS -S /bin/bash
#PBS -q newest 
#PBS -N mat_mul
#PBS -l nodes=1:ppn=10
#PBS -l walltime=00:05:00

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=10
./mat_mul
