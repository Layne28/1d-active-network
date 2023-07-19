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
    in_dir="${scratch_dir}va=${va}/tau=${tau}/lambda=${lambda}/Nx=${Nx}/nx=${nx}/seed=${i}/prod/trajectory"
    ref_file="${scratch_dir}va=${va}/tau=${tau}/lambda=${lambda}/Nx=${Nx}/nx=${nx}/seed=${i}/equil/trajectory/0.dump"
    out_dir="${scratch_dir}va=${va}/tau=${tau}/lambda=${lambda}/Nx=${Nx}/nx=${nx}/seed=${i}/prod/processed_data"

    sbatch -J "get_space_corr_va=${va}_tau=${tau}_lambda=${lambda}_seed=${i}" -o "log/get_space_corr_va=${va}_tau=${tau}_lambda=${lambda}_seed=${i}.o%j" -e "log/get_space_corr_va=${va}_tau=${tau}_lambda=${lambda}_seed=${i}.e%j" scripts/submit_get_space_corr.sh ${in_dir} ${ref_file} ${out_dir} 
done
