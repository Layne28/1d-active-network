#!/bin/bash

scratch_dir="/work/laynefrechette/active-network-results/"

Nx=$1
Ny=$2
Nz=$3
nx=$4
the_type=$5 #disp or vel

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
            the_dir="${scratch_dir}va=${va}/tau=${tau}/lambda=${lambda}/Nx=${Nx}_Ny=${Ny}_Nz=${Nz}/nx=${nx}/"

            sbatch -J "avg_${the_type}_mode_time_corr_va=${va}_tau=${tau}_lambda=${lambda}" -o "log/avg_${the_type}_mode_time_corr_va=${va}_tau=${tau}_lambda=${lambda}.o%j" -e "log/avg_${the_type}_mode_time_corr_va=${va}_tau=${tau}_lambda=${lambda}.e%j" scripts/submit_avg_mode_time_corr.sh ${the_dir} ${the_type}
            #sleep 30
        done
    done
done
