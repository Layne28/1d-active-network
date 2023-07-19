#!/bin/bash

scratch_dir="/work/laynefrechette/passive-network-results/"

Nx=$1
Ny=$2
Nz=$3
kT=$4
nseed=$5
the_type=$6 #disp or vel

for (( i=1; i<=$nseed; i++ ))
do
    the_dir="${scratch_dir}kT=${kT}/Nx=${Nx}_Ny=${Ny}_Nz=${Nz}/seed=${i}/prod/processed_data/"

    sbatch -J "get_${the_type}_mode_time_corr_seed=${i}" -o "log/get_${the_type}_mode_time_corr_seed=${i}.o%j" -e "log/get_${the_type}_mode_time_corr_seed=${i}.e%j" scripts/submit_get_mode_time_corr.sh ${the_dir} ${the_dir} ${the_type}
done
