#!/bin/bash

N=100
nx=10000
taus=(0.030000 0.100000 0.300000 1.000000 10.000000)
lambdas=(0.100000 0.300000 1.000000 3.000000 10.000000)
vas=(0.010000 0.030000 0.100000 0.300000 1.000000)

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for va in "${vas[@]}"
        do
            filename="seeds/seeds_N=${N}_nx=${nx}_va=${va}_tau=${tau}_lambda=${lambda}.txt"
            touch $filename
            for i in {1..5}
            do
                printf "$RANDOM\n" >> $filename
            done
        done
    done
done
