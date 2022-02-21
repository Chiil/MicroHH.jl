#!/bin/bash
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=4
#SBATCH --partition=thin
#SBATCH --time=01:00:00
 
# Load modules for MPI and other parallel libraries
module load 2021
module load foss/2021a

srun julia --project -O3 -t4 drycbl_run.jl --use-mpi --npx 32 --npy 64

