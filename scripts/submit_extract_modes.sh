#!/bin/bash
#Extract Fourier modes of displacement and velocity fields

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH --time=30:00:00
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --exclude=compute-3-17

module load share_modules/ANACONDA/5.3_py3

h5_file=$1
ref_file=$2
out_folder=$3

echo "Running mode extraction in folder: $in_folder"

run_dir="/home/laynefrechette/active-network/scripts/"

python3 $run_dir/extract_modes_h5.py $h5_file $ref_file $out_folder

