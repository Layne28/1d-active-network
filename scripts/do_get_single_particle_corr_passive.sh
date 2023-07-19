#!/bin/bash

scratch_dir="/work/laynefrechette/passive-network-results/"

Nx=$1
Ny=$2
Nz=$3
kT=$4
nseed=$5

for (( i=1; i<=$nseed; i++ ))
do
    h5_file="${scratch_dir}kT=${kT}/Nx=${Nx}_Ny=${Ny}_Nz=${Nz}/seed=${i}/prod/traj.h5"

    ref_file="${scratch_dir}kT=${kT}/Nx=${Nx}_Ny=${Ny}_Nz=${Nz}/seed=${i}/equil/traj.h5"

    out_dir="${scratch_dir}kT=${kT}/Nx=${Nx}_Ny=${Ny}_Nz=${Nz}/seed=${i}/prod/processed_data/"

    sbatch -J "get_${the_type}_single_particle_corr_seed=${i}" -o "log/get_${the_type}_single_particle_corr_seed=${i}.o%j" -e "log/get_${the_type}_single_particle_corr_va=${va}_tau=${tau}_lambda=${lambda}_seed=${i}.e%j" scripts/submit_get_single_particle_corr.sh ${h5_file} ${ref_file} ${out_dir}
done
