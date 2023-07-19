#!/bin/bash
#Average Fourier mode correlations over different trajectories

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH --time=2:00:00
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --exclude=compute-3-17

module load share_modules/ANACONDA/5.3_py3

in_folder=$1
the_type=$2

echo "Running mode correlation average in folder: $in_folder"

run_dir="/home/laynefrechette/active-network/scripts/"

python3 $run_dir/avg_mode_time_corr.py $in_folder $the_type

