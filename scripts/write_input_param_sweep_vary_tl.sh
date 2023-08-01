#!/bin/bash

#N=$1
#nx=$2
#kT=$3
potential=$1
basedir=$2 #work or scratch0
do_subtract_com=$3
do_pin_node=$4

N=100
nx=10000
dt=0.00025
freq=4000
equil_steps=4000000
production_steps=0
kT=0.0

taus=(0.03 0.10 0.30 1.00 10.00)
lambdas=(0.10 0.30 1.00 3.00 10.00)
vas=(0.01 0.03 0.10 0.30 1.00)

com_method='com_unconstrained'
if [ do_subtract_com=='1' ] && [ do_pin_node=='0' ]
then
    com_method='com_fixed'
fi
if [ do_pin_node=='1' ]
then
    com_method='pinned'
fi

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for va in "${vas[@]}"
        do
            python3 scripts/write_input_file.py input_files --va $va --tau $tau --Lambda $lambda --kT $kT --N $N --La $N --nx $nx --network_freq $freq --thermo_freq $freq --dt $dt --equil_steps=$equil_steps --production_steps=$production_steps --potential_type=$potential --do_subtract_com ${do_subtract_com} --do_pin_node ${do_pin_node} --output_dir "/${basedir}/laynefrechette/1d-active-network-results/${com_method}/${potential}/kT=${kT}/va=${va}/tau=${tau}/lambda=${lambda}/N=${N}/nx=${nx}"
        done
    done
done
