#!/bin/bash
#SBATCH --time=0:10:00
#SBATCH --job-name="dbslmm example"
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=5
#SBATCH --output=example.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fredboe@umich.edu  

bash example-dbslmm.sh


