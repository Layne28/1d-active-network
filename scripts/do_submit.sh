#!/bin/bash

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
            input="seeds/seeds_Nx=${Nx}_Ny=${Ny}_Nz=${Nz}_nx=${nx}_va=${va}_tau=${tau}_lambda=${lambda}.txt"
            for (( i=1; i<=$nseed; i++ ))
            do
                echo $i
                sbatch -J "active_net_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_Ny=${Ny}_Nz=${Nz}_nx=${nx}_seed=${i}" -o "log/active_net_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_Ny=${Ny}_Nz=${Nz}_nx=${nx}_seed=${i}.o%j" -e "log/active_net_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_Ny=${Ny}_Nz=${Nz}_nx=${nx}_seed=${i}.e%j" scripts/submit_sims.sh "/home/laynefrechette/active-network/input_files/active_net_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_Ny=${Ny}_Nz=${Nz}_nx=${nx}.in" $i $input
            done
        done
    done
done
