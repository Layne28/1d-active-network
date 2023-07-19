#Compute trajectory average of time autocorrelation functions of displacement field modes

import sys
import os
import glob
import numpy as np
from math import sqrt
import numpy.linalg as la
import h5py

def main():

    in_folder = sys.argv[1]
    out_folder = in_folder + '/avg_data'

    if(os.path.exists(in_folder + '/seed=1/prod/processed_data')):
        data = h5py.File(in_folder + '/seed=1/prod/processed_data/single_particle_t_corr.h5', 'r')
        time = np.array(data['/corr_t/time'])
        a = np.array(data['/parameters/a'])
        edges = np.array(data['/parameters/edges'])
        data.close()
    else:
        print('Input subfolder "seed=1" does not exist. Exiting.')
        exit()        

    if not os.path.exists(out_folder):
        print('Output folder does not exist. Creating it now.')
        os.makedirs(out_folder)

    #Get C(t) for each seed
    cuu_list = []
    cvv_list = []

    for folder in os.listdir(in_folder):
        if 'seed' in folder and folder!='seed_value.txt':
            #print(folder)
            data = h5py.File((in_folder + '/' + folder + '/prod/processed_data/single_particle_t_corr.h5'), 'r')
            cuu_list.append(np.array(data['/corr_t/c_uu']))
            cvv_list.append(np.array(data['/corr_t/c_vv']))

    #Compute trajectory-averaged C(t)
    cuu_avg = np.zeros(np.array(cuu_list[0]).shape)
    cuu_err = np.zeros(np.array(cuu_list[0]).shape)

    cvv_avg = np.zeros(np.array(cvv_list[0]).shape)
    cvv_err = np.zeros(np.array(cvv_list[0]).shape)

    #Get mean C(t)
    for cuu in cuu_list:
        cuu_avg += cuu
    cuu_avg /= len(cuu_list)

    for cvv in cvv_list:
        cvv_avg += cvv
    cvv_avg /= len(cvv_list)

    #Get std error of C(t)
    for cuu in cuu_list:
        cuu_err += (cuu-cuu_avg)**2
    cuu_err = np.sqrt(cuu_err/(len(cuu_list)**2))

    for cvv in cvv_list:
        cvv_err += (cvv-cvv_avg)**2
    cvv_err = np.sqrt(cvv_err/(len(cvv_list)**2))

    #Create output file
    corr_file = h5py.File(out_folder + '/avg_single_particle_t_corr.h5', 'w')
    corr_file.create_dataset('/corr_t/time', data=time)
    corr_file.create_dataset('/parameters/a', data=a)
    corr_file.create_dataset('/parameters/edges', data=edges)
    corr_file.create_dataset('/parameters/ntraj', data=len(cuu_list))

    corr_file.create_dataset('/corr_t/c_uu/corr_avg', data=cuu_avg)
    corr_file.create_dataset('/corr_t/c_uu/corr_err', data=cuu_err)
    corr_file.create_dataset('/corr_t/c_vv/corr_avg', data=cvv_avg)
    corr_file.create_dataset('/corr_t/c_vv/corr_err', data=cvv_err)


main()