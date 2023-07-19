#!/bin/bash

Nx=10
nx=64

#taus=(0.010000 0.030000 0.100000 0.300000 1.000000 10.000000)
#lambdas=(0.010000 0.030000 0.100000 0.300000 1.000000 2.000000)
#vas=(0.010000 0.030000 0.100000 0.300000 1.000000)
#vas=(0.030000 0.100000)
#vas=(0.300000 1.000000)
taus=(1.000000)
lambdas=(1.000000)
vas=(0.100000 0.300000)

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for va in "${vas[@]}"
        do
            sbatch -J "active_net_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_nx=${nx}" -o "active_net_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_nx=${nx}.o%j" -e "active_net_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_nx=${nx}.e%j" scripts/submit_sims.sh "/home/laynefrechette/active-network/conf/active_net_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_nx=${nx}.conf" $RANDOM
        done
    done
done
