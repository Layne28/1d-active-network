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

    #Must be given or assume an equilibrium bond length
    leq = 1.0

    subfolder_list = []
    #Go through all seed subfolders, get data there and append to histogram
    if(os.path.exists(in_folder)):
        for subfolder in sorted(sorted(glob.glob(in_folder + "/seed=*")),key=len):
            print(subfolder)
            subfolder_list.append(subfolder)
    else:
        print('Error: folder does not exist.')
        exit()

    strain_values = []

    for subfolder in subfolder_list:
        #Load data from reference lattice
        eq_traj = h5py.File(subfolder+'/equil/traj.h5')
        N, Lx, Ly, Lz, ref_lat = load_ref(eq_traj)

        #Load trajectory
        traj = h5py.File(subfolder+'/prod/traj.h5')
        pos = np.array(traj['/particles/all/position/value'])
        bonds = np.array(traj['/connectivity/all_bonds'])
        print(pos.shape)
        print(bonds.shape)

        strain_arr = get_strain_arr(pos, bonds, N, Lx, Ly, Lz, leq)

        #Compute strain
        strain_values.append(strain_arr[:,:,-1])

    strain_data = (np.array(strain_values)).flatten()

    #Cumulative distributions, to get percentiles
    strain_sorted = np.sort(strain_data)
    ptiles = np.array(range(strain_data.shape[0]))/float(strain_data.shape[0])

    ptile90_upper = strain_sorted[np.argmax(ptiles>0.9)]
    ptile95_upper = strain_sorted[np.argmax(ptiles>0.95)]
    ptile975_upper = strain_sorted[np.argmax(ptiles>0.975)]
    ptile99_upper = strain_sorted[np.argmax(ptiles>0.99)]

    ptile90_lower = strain_sorted[np.argmax(ptiles>=0.1)]
    ptile95_lower = strain_sorted[np.argmax(ptiles>=0.05)]
    ptile975_lower = strain_sorted[np.argmax(ptiles>=0.025)]
    ptile99_lower = strain_sorted[np.argmax(ptiles>=0.01)]

    with open(in_folder + '/thresh_strain.txt', 'w') as f:
        f.write('percentile strain_upper strain_lower\n')
        f.write('%f %f %f\n' % (0.9, ptile90_upper, ptile90_lower))
        f.write('%f %f %f\n' % (0.95, ptile95_upper, ptile95_lower))
        f.write('%f %f %f\n' % (0.975, ptile975_upper, ptile975_lower))
        f.write('%f %f %f\n' % (0.99, ptile99_upper, ptile99_lower))

    #Histogram
    hist, bin_edges = np.histogram(strain_data, bins=100, density=True)
    bins = (bin_edges[:-1]+bin_edges[1:])/2

    np.savetxt(in_folder + '/strain_hist.txt', np.c_[bins,hist])

    '''
    fig = plt.figure()
    plt.plot(bins, hist)
    plt.axvline(ptile90_upper, linestyle='--', color='black', label='0.90')
    plt.axvline(ptile95_upper, linestyle='--', color='blue', label='0.95')
    plt.axvline(ptile975_upper, linestyle='--', color='purple', label='0.975')
    plt.axvline(ptile99_upper, linestyle='--', color='red', label='0.99')
    plt.axvline(ptile90_lower, linestyle='--', color='black')
    plt.axvline(ptile95_lower, linestyle='--', color='blue')
    plt.axvline(ptile975_lower, linestyle='--', color='purple')
    plt.axvline(ptile99_lower, linestyle='--', color='red')
    plt.legend()
    plt.xlabel(r'$s$')
    plt.xlabel(r'$P(s)$')
    #plt.yscale('log')
    plt.savefig(in_folder + 'strain_histo.png')
    #plt.show()
    '''

def load_ref(traj):

    N = traj['/particles/all/position/value'].shape[1]
    Lx = traj['/particles/all/box/edges'][0]
    Ly = traj['/particles/all/box/edges'][1]
    Lz = traj['/particles/all/box/edges'][2]
    ref_lat = traj['/particles/all/position/value'][0,:]

    return N, Lx, Ly, Lz, ref_lat

#This function takes in a position trajectory and an array of bonds (fixed through trajectory)
@numba.jit(nopython=True)
def get_strain_arr(pos, bond_list, N, Lx, Ly, Lz, leq):

    nframes = pos.shape[0]
    nbonds = bond_list.shape[0]
    strain_arr = np.zeros((nframes, nbonds, 7))

    for t in range(nframes):
        for i in range(nbonds):
            b = bond_list[i,:]
            min_disp_vec = get_min_disp(pos[t,b[0],:], pos[t,b[1],:], Lx, Ly, Lz)
            pos1, pos2 = pos[t,b[0],:], pos[t,b[0],:]-min_disp_vec
            if la.norm(pos1-pos2)>2:
                print('Error!!')
            strain = la.norm(min_disp_vec)-leq
            strain_arr[t][i][0] = pos1[0]
            strain_arr[t][i][1] = pos1[1]
            strain_arr[t][i][2] = pos1[2]
            strain_arr[t][i][3] = pos2[0]
            strain_arr[t][i][4] = pos2[1]
            strain_arr[t][i][5] = pos2[2]
            strain_arr[t][i][6] = strain
    return strain_arr

@numba.jit(nopython=True)
def apply_min_image(disp_r, Lx, Ly, Lz):
    new_disp = np.zeros((disp_r.shape))
    for i in range(new_disp.shape[0]):
        new_disp[i,:] = get_min_disp(disp_r[i,:],np.zeros(3),Lx,Ly,Lz)
    return new_disp

@numba.jit(nopython=True) 
def get_min_disp(r1, r2, Lx, Ly, Lz):
    arr1 = np.array([Lx/2,Ly/2,Lz/2])
    arr2 = np.array([-Lx/2,-Ly/2,-Lz/2])
    arr3 = np.array([Lx, Ly, Lz])
    rdiff = r1-r2
    rdiff = np.where(rdiff>arr1, rdiff-arr3, rdiff)
    rdiff = np.where(rdiff<arr2, rdiff+arr3, rdiff)
    return rdiff

main()