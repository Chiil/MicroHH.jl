#!/bin/bash
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=4
#SBATCH --partition=thin
#SBATCH --time=01:00:00
#SBATCH --exclusive
 
# Load modules for MPI and other parallel libraries
module load 2021
module load foss/2021a

srun julia --project -O3 -t4 drycbl_init.jl --use-mpi --npx 16 --npy 16
srun julia --project -O3 -t4 drycbl_run.jl --use-mpi --npx 16 --npy 16

