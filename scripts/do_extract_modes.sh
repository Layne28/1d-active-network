#!/bin/bash

scratch_dir="/work/laynefrechette/active-network-results/"

Nx=$1
Ny=$2
Nz=$3
nx=$4
nseed=$5

taus=(0.010000 0.100000 0.300000 1.000000 10.000000)
lambdas=(0.010000 0.100000 0.300000 1.000000 2.000000 5.000000)
#vas=(0.001000 0.010000 0.100000 0.300000 1.000000)
vas=(0.100000)

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for va in "${vas[@]}"
        do
            for (( i=1; i<=$nseed; i++ ))
            do
                echo "$tau $lambda $va $i"
                base_dir="${scratch_dir}va=${va}/tau=${tau}/lambda=${lambda}/Nx=${Nx}_Ny=${Ny}_Nz=${Nz}/nx=${nx}/seed=${i}"
                #if [ ! -f "${base_dir}/prod/processed_data/uq.h5" ]; then
                echo "Submitting job."
                mkdir "${base_dir}/prod/processed_data"
                sbatch -J "extract_modes_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_Ny=${Ny}_Nz=${Nz}_nx=${nx}_seed=${i}" -o "log/extract_modes_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_Ny=${Ny}_Nz=${Nz}_nx=${nx}_seed=${i}.o%j" -e "log/extract_modes_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_Ny=${Ny}_Nz=${Nz}_nx=${nx}_seed=${i}.e%j" scripts/submit_extract_modes.sh "${base_dir}/prod/traj.h5" "${base_dir}/equil/traj.h5" "${base_dir}/prod/processed_data"
                #fi
            done
            sleep 30
        done
    done
done
