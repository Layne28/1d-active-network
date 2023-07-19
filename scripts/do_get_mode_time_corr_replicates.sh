#!/bin/bash

scratch_dir="/scratch0/laynefrechette/active-network-results/"
#scratch_dir="/work2/laynefrechette/active-network-results/"

Nx=$1
nx=$2
va=$3
tau=$4
lambda=$5
the_type=$6 #disp or vel

for i in {1..50}
do
    the_dir="${scratch_dir}va=${va}/tau=${tau}/lambda=${lambda}/Nx=${Nx}_Ny=${Nx}_Nz=${Nx}/nx=${nx}/seed=${i}/prod/processed_data"

    sbatch -J "get_${the_type}_mode_time_corr_va=${va}_tau=${tau}_lambda=${lambda}_seed=${i}" -o "log/get_${the_type}_mode_time_corr_va=${va}_tau=${tau}_lambda=${lambda}_seed=${i}.o%j" -e "log/get_${the_type}_mode_time_corr_va=${va}_tau=${tau}_lambda=${lambda}_seed=${i}.e%j" scripts/submit_get_mode_time_corr.sh ${the_dir} ${the_dir} ${the_type}
done
