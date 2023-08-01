#!/bin/bash
#Submit active network simulations

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH --time=96:00:00
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --exclude=compute-3-17

conf_file=$1
seed=$2
seed_file=$3

echo "Running 1d active network with input file: $conf_file"

run_dir="/home/laynefrechette/1d-active-network/bin/"

$run_dir/active_network $conf_file $seed $seed_file

