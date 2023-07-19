#!/bin/bash

scratch_dir="/work/laynefrechette/active-network-results/"

Nx=$1
Ny=$2
Nz=$3
nx=$4

taus=(0.010000 0.100000 0.300000 1.000000 10.000000)
lambdas=(0.010000 0.100000 0.300000 1.000000 2.000000 5.000000)
vas=(0.100000)

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for va in "${vas[@]}"
        do
            echo "$tau $lambda $va"
            base_dir="${scratch_dir}va=${va}/tau=${tau}/lambda=${lambda}/Nx=${Nx}_Ny=${Ny}_Nz=${Nz}/nx=${nx}"
            echo "Submitting job."
            sbatch -J "get_strain_histo_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_Ny=${Ny}_Nz=${Nz}_nx=${nx}" -o "log/get_strain_histo_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_Ny=${Ny}_Nz=${Nz}_nx=${nx}.o%j" -e "log/get_strain_histo_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_Ny=${Ny}_Nz=${Nz}_nx=${nx}.e%j" --wait scripts/submit_get_strain_histo.sh "${base_dir}"
        done
    done
done
