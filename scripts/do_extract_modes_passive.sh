#!/bin/bash

scratch_dir="/work/laynefrechette/passive-network-results/"

Nx=$1
Ny=$2
Nz=$3
kT=$4
nseed=$5

for (( i=1; i<=$nseed; i++ ))
do
    echo "$tau $lambda $va $i"
    base_dir="${scratch_dir}kT=${kT}/Nx=${Nx}_Ny=${Ny}_Nz=${Nz}/seed=${i}"
    #if [ ! -f "${base_dir}/prod/processed_data/uq.h5" ]; then
    echo "Submitting job."
    mkdir "${base_dir}/prod/processed_data"
    sbatch -J "extract_modes_seed=${i}" -o "log/extract_modes_seed=${i}.o%j" -e "log/extract_modes_seed=${i}.e%j" scripts/submit_extract_modes.sh "${base_dir}/prod/traj.h5" "${base_dir}/equil/traj.h5" "${base_dir}/prod/processed_data"
    #fi
done
