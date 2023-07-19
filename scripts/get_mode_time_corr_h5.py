#Compute time autocorrelation functions of displacement field modes

import sys
import os
import glob
import numpy as np
from math import sqrt
import numpy.linalg as la
import h5py
import numba

def main():

    in_folder = sys.argv[1]
    out_folder = sys.argv[2]
    thetype = sys.argv[3] #displacement or velocity

    if thetype=='disp':
        name='u'
        data = h5py.File(in_folder + '/uq.h5','r')
        print('Getting time correlation of displacement modes.')
    elif thetype=='vel':
        name='v'
        data = h5py.File(in_folder + '/vq.h5','r')
        print('Getting time correlation of velocity modes.')
    else:
        print('Error: type of mode not recognized.')
        exit()

    tau_max = 25 #max time (in units of tau) to which to compute correlation function

    if not os.path.exists(out_folder):
        print('Output folder does not exist. Creating it now.')
        os.makedirs(out_folder)

    a = np.array(data['/parameters/a'])
    Lx = data['/parameters/edges'][0]
    Ly = data['/parameters/edges'][1]
    Lz = data['/parameters/edges'][2]
    times = np.array(data['/data/time'])
    indices = np.array(data['/indices'])
    uq_all = np.array(data['/data/value'])
    data.close()

    delta_t = times[1]-times[0]
    frame_diff_max = int(tau_max/delta_t) #max number of data points to which to compute correlation function
    print(frame_diff_max)
    print('Loaded data.')

    #Create output file
    corr_file = h5py.File(out_folder + '/%sq_t_corr.h5' % name, 'w')
    corr_file.create_dataset('/corr_t/time', data=delta_t*np.arange(frame_diff_max))
    corr_file.create_dataset('/corr_q/indices', data=indices)
    corr_file.create_dataset('/parameters/a', data=a)
    corr_file.create_dataset('/parameters/edges', data=np.array([Lx,Ly,Lz]))

    #Compute C(t)
    c_q = np.zeros((indices.shape[0],3,3))
    for i in range(indices.shape[0]):
        n1 = indices[i,0]
        n2 = indices[i,1]
        n3 = indices[i,2]
        print(i, n1, n2 ,n3)

        uqmat = uq_all[:,i,:] #first dimension is time, second is space
        c_t = get_corr_t(uqmat, frame_diff_max)

        corr_file.create_dataset('/corr_t/nx=%d_ny=%d_nz=%d' % (n1, n2, n3), data=c_t) #3x3 matrix for each time
        c_q[i,:,:] = c_t[0,:,:]

    corr_file.create_dataset('/corr_q/corr', data=c_q)

@numba.jit(nopython=True) 
def get_corr_t(uqmat, tmax):
    c_t = np.zeros((tmax,3,3),dtype=np.complex128)
    for t in range(tmax):
        #take outer product of uq(t) and uq(t+delta_t) to get 3x3 matrix, and average over t
        for i in range(3):
            for j in range(3):
                c_t[t][i][j]= np.mean(np.conj(uqmat[:-tmax,i])*uqmat[t:(-tmax+t),j])
        #c_t[t]= np.mean(np.conj(uqmat[:-tmax,:, None])*uqmat[t:(-tmax+t),None,:], axis=0)
    return c_t

main()