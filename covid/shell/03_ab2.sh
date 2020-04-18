#!/bin/bash
#SBATCH -p batch
#SBATCH -t 0-16:00:00
#SBATCH --constraint=[ivb|hsw]
#SBATCH -n 1
#SBATCH --mem=3000
#SBATCH -o out/ab2-0_1.out
module load r/3.6.1-python3
srun Rscript --vanilla main_triton.R 10 1000 3 0.1
