#!/bin/bash

Nx=$1
Ny=$2
Nz=$3
kT=$4

filename="seeds/seeds_Nx=${Nx}_Ny=${Ny}_Nz=${Nz}_kT=${kT}.txt"
touch $filename
for i in {1..100}
do
    printf "$RANDOM\n" >> $filename
done