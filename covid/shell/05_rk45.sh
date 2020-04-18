#!/bin/bash
#SBATCH -p batch
#SBATCH -t 5-00:00:00
#SBATCH --constraint=[ivb|hsw]
#SBATCH -n 1
#SBATCH --mem=3000
#SBATCH -o out/rk45.out
module load r/3.6.1-python3
srun Rscript --vanilla main_triton.R 10 1000 5 0.0
