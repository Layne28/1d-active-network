#Compute time autocorrelation functions of displacement field modes

import sys
import os
import glob
import numpy as np
from math import sqrt
import numpy.linalg as la

def main():

    in_folder = sys.argv[1]
    out_folder = sys.argv[2]
    thetype = sys.argv[3] #displacement or velocity

    if thetype=='disp':
        name='u'
        print('Getting time correlation of displacement modes.')
    elif thetype=='vel':
        name='v'
        print('Getting time correlation of velocity modes.')
    else:
        print('Error: type of mode not recognized.')
        exit()

    n_dt_max = 25000 #max no. of timesteps to which to compute correlation functions (for dt=2e-3 this is 3tau)

    if(os.path.exists(in_folder)):
        mylist = []
        for thefile in sorted(sorted(glob.glob(in_folder + "/%sq*.txt" % name)), key=len):
            print(thefile)
            mylist.append(thefile)
    else:
        print('Input folder does not exist. Exiting.')
        exit()

    if not os.path.exists(out_folder):
        print('Output folder does not exist. Creating it now.')
        os.makedirs(out_folder)

    #Extract value of Nx
    splitlist = in_folder.split('/')
    Ns = [i for i in splitlist if i.startswith('Nx=')]
    Nx = int(((Ns[0].split('_'))[0].split('='))[-1])
    print('Nx: ', Nx)

    #Loop through files
    times = []
    uqs = []
    for thefile in mylist:
        uqs.append(np.loadtxt(thefile,skiprows=1,dtype=complex))
        with open(thefile) as f:
            times.append(int(f.readline().split()[-1]))
        #print(thefile)
    times = np.array(times)
    print(times)
    dt_diff = times[1]-times[0]
    frame_diff_max = n_dt_max//dt_diff #max number of data points to which to compute correlation function
    print(frame_diff_max)
    print('Loaded data.')

    #Get indices of allowed qs
    #assumes lattice constant a=sqrt(2)
    qs = uqs[0][:,:3]
    indices = np.rint(np.real(qs)*(sqrt(2.0)/(2*np.pi))*Nx).astype(int)
    #print(indices)
    for i in range(indices.shape[0]):
        print(qs[i,:], indices[i,:])

    #Compute C(t)
    c_t = {}
    for i in range(indices.shape[0]):
        n1 = indices[i,0]
        n2 = indices[i,1]
        n3 = indices[i,2]
        uqmat = np.array([uq[i,3:] for uq in uqs])
        print(i, uqmat.shape)
        c_t[(n1, n2, n3)] = []
        c_0 = np.mean(np.conj(uqmat[:,:, None])*uqmat[:,None,:], axis=0) #zero-time correlation

        for j in range(frame_diff_max):
            #take outer product of uq(t) and uq(t+delta_t) to get 3x3 matrix, and average over t
            result = np.mean(np.conj(uqmat[:-frame_diff_max,:, None])*uqmat[j:(-frame_diff_max+j),None,:], axis=0)
            c_t[(n1,n2,n3)].append(result)
        np.savez(out_folder + '/%sq_t_corr_nx=%d_ny=%d_nz=%d.npz' % (name, n1, n2, n3), time=dt_diff*np.arange(frame_diff_max), corr=np.array(c_t[(n1,n2,n3)]),corr0=c_0)

main()