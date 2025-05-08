#!/bin/bash
#SBATCH --job-name=parsinv_ex2
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err
#SBATCH --mail-user=esmail.abdulfattah@kaust.edu.sa
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
#SBATCH --time=10:00:00
#SBATCH --report-bindings

#load modules:
module load mpi

#OpenMP settings:
export OMP_NUM_THREADS=1

#run the application:
srun ./main -ni 1000 -lr 1.0 -dr 0.8 -nr 10