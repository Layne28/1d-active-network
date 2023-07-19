#!/bin/bash

scratch_dir="/scratch0/laynefrechette/active-network-results/"
#scratch_dir="/work2/laynefrechette/active-network-results/"

Nx=$1
nx=$2
va=$3
tau=$4
lambda=$5

for i in {1..50}
do
    base_dir="${scratch_dir}va=${va}/tau=${tau}/lambda=${lambda}/Nx=${Nx}_Ny=${Nx}_Nz=${Nx}/nx=${nx}/seed=${i}"
    mkdir "${base_dir}/prod/processed_data"
    sbatch -J "extract_modes_va=${va}_tau=${tau}_lambda=${lambda}_seed=${i}" -o "log/extract_modes_va=${va}_tau=${tau}_lambda=${lambda}_seed=${i}.o%j" -e "log/extract_modes_va=${va}_tau=${tau}_lambda=${lambda}_seed=${i}.e%j" scripts/submit_extract_modes.sh "${base_dir}/equil/trajectory/0.dump" "${base_dir}/prod/trajectory" "${base_dir}/prod/processed_data"
done
