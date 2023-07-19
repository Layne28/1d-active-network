#!/bin/bash

scratch_dir="/scratch0/laynefrechette/active-network-results/"

Nx=$1
nx=$2

taus=(0.010000 0.100000 1.000000 10.000000)
lambdas=(0.010000 0.100000 1.000000 2.000000 10.000000)
vas=(0.001000 0.010000 0.100000 0.300000 1.000000)

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for va in "${vas[@]}"
        do
            the_dir="${scratch_dir}va=${va}/tau=${tau}/lambda=${lambda}/Nx=${Nx}/nx=${nx}/prod/processed_data"

            sbatch -J "get_mode_time_corr_va=${va}_tau=${tau}_lambda=${lambda}" -o "get_mode_time_corr_va=${va}_tau=${tau}_lambda=${lambda}.o%j" -e "get_mode_time_corr_va=${va}_tau=${tau}_lambda=${lambda}.e%j" scripts/submit_get_mode_time_corr.sh ${the_dir} ${the_dir}
        done
    done
done
