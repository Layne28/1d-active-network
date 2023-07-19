#!/bin/bash

Nx=10
nx=64

taus=(0.010000 0.030000 0.100000 0.300000 1.000000 10.000000)
lambdas=(0.010000 0.030000 0.100000 0.300000 1.000000 2.000000)
vas=(0.010000 0.030000 0.100000 0.300000 1.000000)

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for va in "${vas[@]}"
        do
            python3 scripts/write_conf_file.py conf --va $va --tau $tau --Lambda $lambda --Nx $Nx --Ny $Nx --Nz $Nx --nx $nx --ny $nx --nz $nx --output_dir "/scratch0/laynefrechette/active-network-results/va=${va}/tau=${tau}/lambda=${lambda}/Nx=${Nx}"
        done
    done
done
