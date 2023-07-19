#!/bin/bash

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH --time=00:10:00
#SBATCH -N 1 
#SBATCH -n 1

module load share_modules/ANACONDA/5.3_py3

python3 scripts/test.py $1