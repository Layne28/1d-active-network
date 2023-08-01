#!/bin/bash

com_method=$1

N=100
nx=10000
nseed=1
potential_type='fene'

taus=(0.030000 0.100000 0.300000 1.000000 10.000000)
lambdas=(0.100000 0.300000 1.000000 3.000000 10.000000)
vas=(0.010000 0.030000 0.100000 0.300000 1.000000)

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for va in "${vas[@]}"
        do
            input="seeds/seeds_N=${N}_nx=${nx}_va=${va}_tau=${tau}_lambda=${lambda}.txt"
            for (( i=1; i<=$nseed; i++ ))
            do
                echo $i
                sbatch -J "${potential_type}_${com_method}_net_va=${va}_tau=${tau}_lambda=${lambda}_N=${N}_nx=${nx}_seed=${i}" -o "log/${potential_type}_${com_method}_net_va=${va}_tau=${tau}_lambda=${lambda}_N=${N}_nx=${nx}_seed=${i}.o%j" -e "log/${potential_type}_${com_method}_net_va=${va}_tau=${tau}_lambda=${lambda}_N=${N}_nx=${nx}_seed=${i}.e%j" scripts/submit_sims.sh "/home/laynefrechette/1d-active-network/input_files/${com_method}/${potential_type}_net_va=${va}_tau=${tau}_lambda=${lambda}_N=${N}_nx=${nx}.in" $i $input
            done
        done
    done
done
