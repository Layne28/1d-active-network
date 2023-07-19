#!/bin/bash

scratch_dir="/scratch0/laynefrechette/active-network-results/"

Nx=$1
nx=$2
va=$3
tau=$4
lambda=$5
i=$6

the_dir="${scratch_dir}va=${va}/tau=${tau}/lambda=${lambda}/Nx=${Nx}/nx=${nx}/seed=${i}/prod/processed_data"

sbatch -J "get_mode_traj_va=${va}_tau=${tau}_lambda=${lambda}_seed=${i}" -o "get_mode_traj_va=${va}_tau=${tau}_lambda=${lambda}_seed=${i}.o%j" -e "get_mode_traj_va=${va}_tau=${tau}_lambda=${lambda}_seed=${i}.e%j" scripts/submit_get_mode_traj.sh ${the_dir} ${the_dir}
