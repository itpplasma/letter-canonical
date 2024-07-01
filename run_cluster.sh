#!/bin/bash
#SBATCH --job-name=letter-canonical
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --output=output_%j.txt

export OMP_NUM_THREADS=64
export MKL_NUM_THREADS=1

./letter-canonical.x
