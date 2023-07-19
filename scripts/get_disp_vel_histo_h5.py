#Get displacement histogram from collection of h5 files

import numpy as np
import numpy.linalg as la
import pylab as plt
import sys
import os
import glob
import h5py

def main():

    in_folder = sys.argv[1]

    subfolder_list = []
    #Go through all seed subfolders, get data there and append to histogram
    
    if(os.path.exists(in_folder)):
        for subfolder in sorted(sorted(glob.glob(in_folder + "/seed=*")),key=len):
            print(subfolder)
            subfolder_list.append(subfolder)
    else:
        print('Error: folder does not exist.')
        exit()
    
    for subfolder in sorted(sorted(glob.glob(in_folder + "/seed=*")),key=len):
        print(subfolder)
        subfolder_list.append(subfolder)
    
    disp_values = []
    vel_values = []

    for subfolder in subfolder_list:
        #Load data from reference lattice
        eq_traj = h5py.File(subfolder+'/equil/traj.h5', 'r')
        N, Lx, Ly, Lz, ref_lat = load_ref(eq_traj)
        eq_traj.close()

        #Load trajectory
        traj = h5py.File(subfolder+'/prod/traj.h5', 'r')
        disp_values.append(np.array(traj['/particles/all/position/value'])-ref_lat)
        vel_values.append(np.array(traj['/particles/all/velocity/value']))
        traj.close()

    disp_values = np.array(disp_values)
    vel_values = np.array(vel_values)
    
    ux_data = (disp_values[:,:,:,0]).flatten()
    uy_data = (disp_values[:,:,:,1]).flatten()
    uz_data = (disp_values[:,:,:,2]).flatten()
    mag_data = (la.norm(disp_values,axis=3)).flatten()

    vx_data = (vel_values[:,:,:,0]).flatten()
    vy_data = (vel_values[:,:,:,1]).flatten()
    vz_data = (vel_values[:,:,:,2]).flatten()
    speed_data = (la.norm(vel_values,axis=3)).flatten()

    #Cumulative distributions, to get percentiles
    mag_sorted = np.sort(mag_data)
    speed_sorted = np.sort(speed_data)
    ptiles = np.array(range(mag_data.shape[0]))/float(mag_data.shape[0])

    ptile90_mag = mag_sorted[np.argmax(ptiles>0.9)]
    ptile95_mag = mag_sorted[np.argmax(ptiles>0.95)]
    ptile99_mag = mag_sorted[np.argmax(ptiles>0.99)]

    ptile90_speed = speed_sorted[np.argmax(ptiles>0.9)]
    ptile95_speed = speed_sorted[np.argmax(ptiles>0.95)]
    ptile99_speed = speed_sorted[np.argmax(ptiles>0.99)]

    with open(in_folder + '/thresh_disp_vel.txt', 'w') as f:
        f.write('percentile disp vel\n')
        f.write('%f %f %f\n' % (0.9, ptile90_mag, ptile90_speed))
        f.write('%f %f %f\n' % (0.95, ptile95_mag, ptile95_speed))
        f.write('%f %f %f\n' % (0.99, ptile99_mag, ptile99_speed))

    #Histograms
    hist_ux, bin_edges_ux = np.histogram(ux_data, bins=100, density=True)
    bins_ux = (bin_edges_ux[:-1]+bin_edges_ux[1:])/2
    hist_uy, bin_edges_uy = np.histogram(uy_data, bins=100, density=True)
    bins_uy = (bin_edges_uy[:-1]+bin_edges_uy[1:])/2
    hist_uz, bin_edges_uz = np.histogram(uz_data, bins=100, density=True)
    bins_uz = (bin_edges_uz[:-1]+bin_edges_uz[1:])/2

    hist_mag, bin_edges_mag = np.histogram(mag_data, bins=100, density=True)
    bins_mag = (bin_edges_mag[:-1]+bin_edges_mag[1:])/2

    hist_vx, bin_edges_vx = np.histogram(vx_data, bins=100, density=True)
    bins_vx = (bin_edges_vx[:-1]+bin_edges_vx[1:])/2
    hist_vy, bin_edges_vy = np.histogram(vy_data, bins=100, density=True)
    bins_vy = (bin_edges_vy[:-1]+bin_edges_vy[1:])/2
    hist_vz, bin_edges_vz = np.histogram(vz_data, bins=100, density=True)
    bins_vz = (bin_edges_vz[:-1]+bin_edges_vz[1:])/2

    hist_speed, bin_edges_speed = np.histogram(speed_data, bins=100, density=True)
    bins_speed = (bin_edges_speed[:-1]+bin_edges_speed[1:])/2

    np.savetxt(in_folder + '/ux_hist.txt', np.c_[bins_ux,hist_ux])
    np.savetxt(in_folder + '/uy_hist.txt', np.c_[bins_uy,hist_uy])
    np.savetxt(in_folder + '/uz_hist.txt', np.c_[bins_uz,hist_uz])
    np.savetxt(in_folder + '/umag_hist.txt', np.c_[bins_mag,hist_mag])

    np.savetxt(in_folder + '/vx_hist.txt', np.c_[bins_vx,hist_vx])
    np.savetxt(in_folder + '/vy_hist.txt', np.c_[bins_vy,hist_vy])
    np.savetxt(in_folder + '/vz_hist.txt', np.c_[bins_vz,hist_vz])
    np.savetxt(in_folder + '/speed_hist.txt', np.c_[bins_speed,hist_speed])

    '''
    fig = plt.figure()
    plt.plot(bins_ux, hist_ux, label='ux')
    plt.plot(bins_uy, hist_uy, label='uy')
    plt.plot(bins_uz, hist_uz, label='uz')
    #plt.yscale('log')
    plt.savefig(in_folder + 'disp_component_histo.png')

    fig = plt.figure()
    plt.plot(bins_mag, hist_mag)
    plt.axvline(ptile90_mag, linestyle='--', color='black', label='0.90')
    plt.axvline(ptile95_mag, linestyle='--', color='blue', label='0.95')
    plt.axvline(ptile99_mag, linestyle='--', color='red', label='0.99')
    plt.legend()
    #plt.yscale('log')
    plt.savefig(in_folder + 'disp_mag_histo.png')

    fig = plt.figure()
    plt.plot(bins_vx, hist_vx, label='vx')
    plt.plot(bins_vy, hist_vy, label='vy')
    plt.plot(bins_vz, hist_vz, label='vz')
    plt.savefig(in_folder + 'vel_component_histo.png')

    fig = plt.figure()
    plt.plot(bins_speed, hist_speed)
    plt.axvline(ptile90_speed, linestyle='--', color='black', label='0.90')
    plt.axvline(ptile95_speed, linestyle='--', color='blue', label='0.95')
    plt.axvline(ptile99_speed, linestyle='--', color='red', label='0.99')
    plt.legend()
    plt.savefig(in_folder + 'speed_histo.png')
    '''

    #plt.show()

def load_ref(traj):

    N = traj['/particles/all/position/value'].shape[1]
    Lx = traj['/particles/all/box/edges'][0]
    Ly = traj['/particles/all/box/edges'][1]
    Lz = traj['/particles/all/box/edges'][2]
    ref_lat = np.array(traj['/particles/all/position/value'][0,:])

    return N, Lx, Ly, Lz, ref_lat

main()