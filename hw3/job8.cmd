#PBS -S /bin/bash
#PBS -q newest
#PBS -N job8
#PBS -l nodes=1:ppn=8
#PBS -l walltime=00:05:00

cd $PBS_O_WORKDIR 
export OMP_NUM_THREADS=8
./hello_omp_grape
