#!/bin/bash
#Get spatial correlations of displacements and velocities

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH --time=10:00:00
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --exclude=compute-3-17

module load share_modules/ANACONDA/5.3_py3

in_folder=$1
ref_file=$2
out_folder=$3

echo "Running space correlation in folder: $in_folder"

run_dir="/home/laynefrechette/active-network/scripts/"

python3 $run_dir/get_space_corr.py $in_folder $ref_file $out_folder

