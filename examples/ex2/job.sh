#!/bin/bash
#SBATCH --job-name=parsinv_ex2
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err
#SBATCH --mail-user=djidenou.montcho@kaust.edu.sa
#SBATCH --mail-type=ALL
#SBATCH --partition=workq
#SBATCH --constraint=intel
#SBATCH --nodes=10
#SBATCH --ntasks=200
#SBATCH --ntasks-per-node=20
#SBATCH --ntasks-per-socket=10
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --cpu-bind=threads
#SBATCH --hint=nomultithread
#SBATCH --mem-bind=v,local
#SBATCH --mem-per-cpu=2G
#SBATCH --time=15:00:00
#SBATCH --report-bindings

#load modules:
module load openmpi

#OpenMP settings:
export OMP_NUM_THREADS=1
export PETSC_DIR=/home/montchd/test_abylay/petsc
export PETSC_ARCH=arch-linux-c-debug
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PETSC_DIR/$PETSC_ARCH/lib

#run the application:
srun ./main -ni 100 -lr 1.0 -dr 0.8 -nr 10
