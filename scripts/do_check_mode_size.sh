#!/bin/bash

scratch_dir="/scratch0/laynefrechette/active-network-results/"

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
                base_dir="${scratch_dir}va=${va}/tau=${tau}/lambda=${lambda}/Nx=${Nx}_Ny=${Ny}_Nz=${Nz}/nx=${nx}/seed=${i}/prod/processed_data"
                python3 scripts/check_mode_size.py "${base_dir}"
            done
        done
    done
done
