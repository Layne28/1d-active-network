#!/bin/bash
#Get strain histogram

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH --time=10:00:00
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --exclude=compute-3-17

module load share_modules/ANACONDA/5.3_py3

in_folder=$1

echo "Running strain histogram calculation in folder: $in_folder"

run_dir="/home/laynefrechette/active-network/scripts/"

python3 $run_dir/get_strain_histo_h5.py $in_folder

