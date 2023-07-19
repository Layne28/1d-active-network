#!/bin/bash
#Get time correlation of Fourier modes

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH --time=10:00:00
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --exclude=compute-3-17

module load share_modules/ANACONDA/5.3_py3

in_folder=$1
out_folder=$2
the_type=$3

echo "Running mode time correlation in folder: $in_folder"

run_dir="/home/laynefrechette/active-network/scripts/"

python3 $run_dir/get_mode_time_corr_h5.py $in_folder $out_folder $the_type

