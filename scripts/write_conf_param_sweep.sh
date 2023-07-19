#!/bin/bash

Nx=$1
Ny=$2
Nz=$3
nx=$4

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
            python3 scripts/write_conf_file.py conf --va $va --tau $tau --Lambda $lambda --Nx $Nx --Ny $Ny --Nz $Nz --Lax $Nx --Lay $Ny --Laz $Nz --nx $nx --ny $nx --nz $nx --network_freq 100 --thermo_freq 100 --output_dir "/work/laynefrechette/active-network-results/va=${va}/tau=${tau}/lambda=${lambda}/Nx=${Nx}_Ny=${Ny}_Nz=${Nz}/nx=${nx}"
        done
    done
done
