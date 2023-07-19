#Compute "single particle" (r=0) time autocorrelation functions of displacements and velocities

import sys
import os
import glob
import numpy as np
from math import sqrt
import numpy.linalg as la
import h5py
import numba

def main():

    h5_file = sys.argv[1] #trajectory file in h5md format
    ref_file = sys.argv[2] #trajectory file containing reference lattice
    out_folder = sys.argv[3] #where to put correlation functions

    #Load data from reference lattice
    eq_traj = h5py.File(ref_file)
    N, Lx, Ly, Lz, a, ref_lat = load_ref(eq_traj)
    #print(ref_lat)
    eq_traj.close()
    print('loaded reference configuration')

    #Load trajectory
    traj = h5py.File(h5_file)
    pos_all = np.array(traj['/particles/all/position/value'])
    pos_all_unwrapped = get_pos_unwrapped(pos_all, np.array(traj['/particles/all/image/value']), Lx, Ly, Lz)
    vel_all = np.array(traj['/particles/all/velocity/value'])
    tsteps = np.array(traj['/particles/all/position/step'])
    times = np.array(traj['/particles/all/position/time'])
    traj.close()

    print('loaded trajectory')

    disp_all = pos_all_unwrapped - ref_lat

    tau_max = 20 #max time (in units of tau) to which to compute correlation function

    if not os.path.exists(out_folder):
        print('Output folder does not exist. Creating it now.')
        os.makedirs(out_folder)

    delta_t = times[1]-times[0]
    frame_diff_max = int(tau_max/delta_t) #max number of data points to which to compute correlation function
    print(frame_diff_max)
    print('Loaded data.')

    #Create output file
    corr_file = h5py.File(out_folder + '/single_particle_t_corr.h5', 'w')
    corr_file.create_dataset('/corr_t/time', data=delta_t*np.arange(frame_diff_max))
    corr_file.create_dataset('/parameters/a', data=a)
    corr_file.create_dataset('/parameters/edges', data=np.array([Lx,Ly,Lz]))

    #Compute time correlation function
    c_uu = get_corr_t(disp_all, frame_diff_max)
    c_vv = get_corr_t(vel_all, frame_diff_max)

    corr_file.create_dataset('/corr_t/c_uu', data=c_uu) #3x3 matrix for each time
    corr_file.create_dataset('/corr_t/c_vv', data=c_vv) #3x3 matrix for each time

@numba.jit(nopython=True) 
def get_corr_t(data, tmax):
    c_t = np.zeros((tmax,3,3))
    for t in range(tmax):
        #take outer product of data(t) and data(t+delta_t) to get 3x3 matrix, and average over t
        for i in range(3):
            for j in range(3):
                c_t[t][i][j]= np.mean(np.conj(data[:-tmax,i])*data[t:(-tmax+t),j])
    return c_t

@numba.jit(nopython=True)
def get_pos_unwrapped(pos, image, Lx, Ly, Lz):

    shift = np.zeros(image.shape)
    shift[:,:,0] = image[:,:,0]*Lx
    shift[:,:,1] = image[:,:,1]*Ly
    shift[:,:,2] = image[:,:,2]*Lz

    pos_unwrapped = pos + shift
    
    return pos_unwrapped

def load_ref(traj):

    N = int(traj['/particles/all/position/value'].shape[1])
    Lx = float(traj['/particles/all/box/edges'][0])
    Ly = float(traj['/particles/all/box/edges'][1])
    Lz = float(traj['/particles/all/box/edges'][2])
    a = traj['/particles/all/box/lattice_constant']
    a = np.array(a)
    #print(a.dtype)
    ref_lat = np.array(traj['/particles/all/position/value'][0,:])

    return N, Lx, Ly, Lz, a, ref_lat

main()