#Compute trajectory average of time autocorrelation functions of displacement field modes

import sys
import os
import glob
import numpy as np
from math import sqrt
import numpy.linalg as la

def main():

    in_folder = sys.argv[1]
    out_folder = in_folder + '/avg_data'
    thetype = sys.argv[2] #displacement or velocity

    if thetype=='disp':
        name='u'
        print('Getting spatial correlation of displacements.')
    elif thetype=='vel':
        name='v'
        print('Getting spatial correlation of velocities.')

    if not os.path.exists(in_folder):
        print('Input folder does not exist. Exiting.')
        exit()

    if not os.path.exists(out_folder):
        print('Output folder does not exist. Creating it now.')
        os.makedirs(out_folder)

    #Compute trajectory-averaged C(r)
    data_list = []
    for folder in os.listdir(in_folder):
        if 'seed' in folder and folder!='seed_value.txt':
            #print(folder)
            data = np.load((in_folder + '/' + folder + '/prod/processed_data/%s_r_corr.npz') % name)
            data_list.append(data)
    print(len(data_list))
        
    c_r = np.zeros((data_list[0]['pos'].shape[0],3,3),dtype='complex')
    c_r_err = np.zeros((data_list[0]['pos'].shape[0],3,3),dtype='complex')

    #Get mean C(r)
    for data in data_list:
        c_r += data['corr']
    c_r /= len(data_list)
        
    #Get std error of C(r)
    for data in data_list:
        c_r_err += (np.real(data['corr'])-np.real(c_r))**2
    c_r_err = np.sqrt(c_r_err/(len(data_list)**2)) #n_traj^2 bc it's standard error

    np.savez(out_folder + '/%s_r_corr_avg.npz' % name, pos=data_list[0]['pos'], corr=np.array(c_r), stderr=np.array(c_r_err), nsample=len(data_list))
    
main()