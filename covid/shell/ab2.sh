#!/bin/bash
#SBATCH -p batch
#SBATCH -t 0-02:00:00
#SBATCH --constraint=[ivb|hsw]
#SBATCH -n 1
#SBATCH --mem=2500
#SBATCH -o out/ab.out
module load r/3.6.1-python3
srun Rscript --vanilla main_triton.R 10 200 3 0.5
