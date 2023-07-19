#!/bin/bash

dat_dir="/scratch0/laynefrechette/active-network-results/"

taus=(0.010000 0.030000 0.100000 0.300000 1.000000 10.000000)
lambdas=(0.010000 0.030000 0.100000 0.300000 1.000000 2.000000)
vas=(0.030000 0.100000 0.300000)

for tau in "${taus[@]}"
do
    for lambda in "${lambdas[@]}"
    do
        for va in "${vas[@]}"
        do
            prod_dir="${dat_dir}va=${va}/tau=${tau}/lambda=${lambda}/prod"
            equil_dir="${dat_dir}va=${va}/tau=${tau}/lambda=${lambda}/equil"
            if [ ! -d ${equil_dir} ]; then
                echo "Equil dir does not exist. Creating it and moving data there now."
                mkdir ${equil_dir}
                mkdir "${equil_dir}/trajectory"
                mkdir "${equil_dir}/topology"
            fi
            equil_traj="${equil_dir}/trajectory"
            if [ -z "$(ls -A ${equil_traj})" ]; then
                echo "Equil dir is empty. Moving files from prod dir..."
                cp "${prod_dir}/topology/0.dump" "${equil_dir}/topology/"
                mv "${prod_dir}/trajectory"/?.dump "${equil_dir}/trajectory/"
                mv "${prod_dir}/trajectory"/??.dump "${equil_dir}/trajectory/"
                mv "${prod_dir}/trajectory"/???.dump "${equil_dir}/trajectory/"
                mv "${prod_dir}/trajectory"/1???.dump "${equil_dir}/trajectory/"
                mv "${prod_dir}/trajectory"/2???.dump "${equil_dir}/trajectory/"
                mv "${prod_dir}/trajectory"/3???.dump "${equil_dir}/trajectory/"
                mv "${prod_dir}/trajectory"/4???.dump "${equil_dir}/trajectory/"
                echo "Done."
            fi
        done
    done
done
