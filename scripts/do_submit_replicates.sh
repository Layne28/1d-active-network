#!/bin/bash

Nx=$1
Ny=$2
Nz=$3
nx=$4
va=$5
tau=$6
lambda=$7

input="seeds/seeds_Nx=${Nx}_Ny=${Ny}_Nz=${Nz}_nx=${nx}_va=${va}_tau=${tau}_lambda=${lambda}.txt"

for i in {1..5}
do
    sbatch -J "active_net_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_nx=${nx}_seed=${i}" -o "log/active_net_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_nx=${nx}_seed=${i}.o%j" -e "log/active_net_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_nx=${nx}_seed=${i}.e%j" scripts/submit_sims.sh "/home/laynefrechette/active-network/conf/active_net_va=${va}_tau=${tau}_lambda=${lambda}_Nx=${Nx}_Ny=${Ny}_Nz=${Nz}_nx=${nx}.conf" $i $input
done
