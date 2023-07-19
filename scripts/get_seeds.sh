#!/bin/bash

Nx=$1
Ny=$2
Nz=$3
nx=$4

taus=(0.010000 0.100000 0.300000 1.000000 10.000000)
lambdas=(0.010000 0.100000 0.300000 1.000000 2.000000 5.000000)
vas=(0.100000)

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for va in "${vas[@]}"
        do
            #filename="seeds/seeds_Nx=${Nx}_nx=${nx}_va=${va}_tau=${tau}_lambda=${lambda}.txt"
            filename="seeds/seeds_Nx=${Nx}_Ny=${Ny}_Nz=${Nz}_nx=${nx}_va=${va}_tau=${tau}_lambda=${lambda}.txt"
            touch $filename
            for i in {1..100}
            do
                printf "$RANDOM\n" >> $filename
            done
        done
    done
done
