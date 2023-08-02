#Get strain histogram from collection of h5 files

import numba
import numpy as np
import numpy.linalg as la
import pylab as plt
import sys
import os
import glob
import h5py

def main():

    in_folder = sys.argv[1]

    run_tests() #check that strain is computed correctly for a couple special cases

    which_subdir = 'prod' #load from 'equil' or 'prod'
    nchunks=10

    #Must be given or assume an equilibrium bond length
    leq = 1.0

    #Go through all seed subfolders, get data there and append to histogram
    if(not(os.path.exists(in_folder))):
        print('Error: folder does not exist.')
        exit()

    strain_values = []

    #Load data from reference lattice
    eq_traj = h5py.File(in_folder+'/equil/traj.h5')
    N, Lx, ref_lat = load_ref(eq_traj)

    #Load trajectory
    traj = h5py.File(in_folder+'/%s/traj.h5' % which_subdir)
    pos = np.array(traj['/particles/all/position/value'])
    print(pos.shape)

    #Compute strain
    strain_arr = get_strain_arr(pos, N, Lx, leq)

    #Chunk the strain trajectory into segments
    strain_dict = {}
    strain_dict['nchunks'] = nchunks
    seglen = strain_arr.shape[0]//nchunks
    for n in range(nchunks):
        strain_values.append(strain_arr[(n*seglen):((n+1)*seglen),:,-1])

        #Make histograms
        hist, bin_edges = np.histogram(strain_values[n], bins=100, density=True)
        bins = (bin_edges[:-1]+bin_edges[1:])/2

        strain_dict['bins_%d' % n] = bins
        strain_dict['hist_%d' % n] = hist

    np.savez(in_folder + '/strain_hist.npz', **strain_dict)
    container = np.load(in_folder + '/strain_hist.npz')
    e = {name: container[name] for name in container}
    print(e)

def load_ref(traj):

    N = traj['/particles/all/position/value'].shape[1]
    Lx = traj['/particles/all/box/edges'][0]
    ref_lat = traj['/particles/all/position/value'][0,:]

    return N, Lx, ref_lat

@numba.jit(nopython=True)
def get_strain_arr(pos, N, Lx, leq):

    nframes = pos.shape[0]
    nbonds = N
    strain_arr = np.zeros((nframes, nbonds, 3))

    for t in range(nframes):
        for i in range(nbonds):
            min_disp_vec = get_min_disp(pos[t,i,:], pos[t,(i+1)%N,:], Lx)
            pos1, pos2 = pos[t,i,:], pos[t,i,:]-min_disp_vec
            if la.norm(pos1-pos2)>2:
                print('Error!!')
            strain = la.norm(min_disp_vec)-leq
            strain_arr[t][i][0] = pos1[0]
            strain_arr[t][i][1] = pos2[0]
            strain_arr[t][i][2] = strain
    return strain_arr

@numba.jit(nopython=True)
def apply_min_image(disp_r, Lx):
    new_disp = np.zeros((disp_r.shape))
    for i in range(new_disp.shape[0]):
        new_disp[i,:] = get_min_disp(disp_r[i,:],np.zeros(1),Lx)
    return new_disp

@numba.jit(nopython=True) 
def get_min_disp(r1, r2, Lx):
    arr1 = np.array([Lx/2])
    arr2 = np.array([-Lx/2])
    arr3 = np.array([Lx])
    rdiff = r1-r2
    rdiff = np.where(rdiff>arr1, rdiff-arr3, rdiff)
    rdiff = np.where(rdiff<arr2, rdiff+arr3, rdiff)
    return rdiff

def run_tests():
    return 0

main()