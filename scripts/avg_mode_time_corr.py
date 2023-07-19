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
    thetype = sys.argv[2] #displacement or velocity

    if thetype=='disp':
        name='u'
        print('Getting time correlation of displacement modes.')
    elif thetype=='vel':
        name='v'
        print('Getting time correlation of velocity modes.')

    #Get allowed wavevector indices, box parameters, and time
    if(os.path.exists(in_folder + '/seed=1/prod/processed_data')):
        traj = h5py.File(in_folder + '/seed=1/prod/processed_data/%sq_t_corr.h5' % name, 'r')
        indices = np.array(traj['/corr_q/indices'])
        times = np.array(traj['/corr_t/time'])
        a = np.array(traj['/parameters/a'])
        Lx = traj['/parameters/edges'][0]
        Ly = traj['/parameters/edges'][1]
        Lz = traj['/parameters/edges'][2]
        traj.close()
    else:
        print('Input subfolder "seed=1" does not exist. Exiting.')
        exit()        

    if not os.path.exists(out_folder):
        print('Output folder does not exist. Creating it now.')
        os.makedirs(out_folder)

    #Get C(q,t) and C(q) for each seed
    cq_list = []
    ct_list = []

    for folder in os.listdir(in_folder):
        if 'seed' in folder and folder!='seed_value.txt':
            #print(folder)
            data = h5py.File((in_folder + '/' + folder + '/prod/processed_data/%sq_t_corr.h5') % name, 'r')
            cq_list.append(data['/corr_q/corr'])
            ct_list.append(data['/corr_t'])

    #Compute trajectory-averaged C(q)
    cq_avg = np.zeros(np.array(cq_list[0]).shape)
    cq_err = np.zeros(np.array(cq_list[0]).shape)

    #Get mean C(q)
    for cq in cq_list:
        cq_avg += cq
    cq_avg /= len(cq_list)

    #Get std error of C(q)
    for cq in cq_list:
        cq_err += (np.real(cq)-np.real(cq_avg))**2
    cq_err = np.sqrt(cq_err/(len(cq_list)**2))
    
    #Compute trajectory-averaged C(t) for each wavenumber
    ct_avg = {}
    ct_err = {}
    for i in range(indices.shape[0]):
        nx = int(indices[i][0])
        ny = int(indices[i][1])
        nz = int(indices[i][2])
        triple = 'nx=%d_ny=%d_nz=%d' % (nx, ny, nz)
        
        ct_avg[(nx, ny, nz)] = np.zeros((times.shape[0],3,3),dtype='complex')
        ct_err[(nx, ny, nz)] = np.zeros((times.shape[0],3,3),dtype='complex')

        #Get mean C(t)
        for ct in ct_list:
            ct_avg[(nx, ny, nz)] += ct['nx=%d_ny=%d_nz=%d' % (nx, ny, nz)]
        ct_avg[(nx, ny, nz)] /= len(ct_list)
        
        #Get std error of C(t)
        for ct in ct_list:
            #Get real and imaginary contributions to variance
            ct_err[(nx, ny, nz)] += (np.real(ct['nx=%d_ny=%d_nz=%d' % (nx, ny, nz)])-np.real(ct_avg[(nx, ny, nz)]))**2
        ct_err[(nx, ny, nz)] = np.sqrt(ct_err[(nx, ny, nz)]/(len(ct_list)**2)) #n_traj^2 bc it's standard error

    #Create output file
    corr_file = h5py.File(out_folder + '/%sq_t_corr.h5' % name, 'w')
    corr_file.create_dataset('/corr_t/time', data=times)
    corr_file.create_dataset('/corr_q/indices', data=indices)
    corr_file.create_dataset('/parameters/a', data=a)
    corr_file.create_dataset('/parameters/edges', data=np.array([Lx,Ly,Lz]))
    corr_file.create_dataset('/parameters/ntraj', data=len(cq_list))

    corr_file.create_dataset('/corr_q/corr_avg', data=cq_avg)
    corr_file.create_dataset('/corr_q/corr_err', data=cq_err)

    for i in range(indices.shape[0]):
        nx = int(indices[i][0])
        ny = int(indices[i][1])
        nz = int(indices[i][2])
        triple = 'nx=%d_ny=%d_nz=%d' % (nx, ny, nz)
        corr_file.create_dataset('/corr_t/%s/corr_avg' % triple, data=ct_avg[(nx, ny, nz)])
        corr_file.create_dataset('/corr_t/%s/corr_err' % triple, data=ct_err[(nx, ny, nz)])
    
main()