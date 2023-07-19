#Cluster displacement and velocity fields

import numpy as np
import numpy.linalg as la
import pylab as plt
import sys
import os
import glob
import h5py
import pandas as pd
import numba

def main():

    in_folder = sys.argv[1]
    percentile_thresh = float(sys.argv[2])

    run_tests(in_folder)

    subfolder_list = []
    
    if(os.path.exists(in_folder)):
        for subfolder in sorted(sorted(glob.glob(in_folder + "/seed=*")),key=len):
            print(subfolder)
            subfolder_list.append(subfolder)
    else:
        print('Error: folder does not exist.')
        exit()

    #Load thresholds
    thresh_table = pd.read_csv(in_folder + '/thresh_disp_vel.txt', sep=' ')
    print(thresh_table)

    disp_thresh = thresh_table.loc[thresh_table['percentile']==percentile_thresh, 'disp'].item()
    vel_thresh = thresh_table.loc[thresh_table['percentile']==percentile_thresh, 'vel'].item()
    print(disp_thresh, vel_thresh)

    #Initialize data containers
    disp_cluster_sizes = [] #number of nodes in clusters
    disp_cluster_areas = [] #number of surface-exposed nodes in clusters
    all_num_disp_clusters = np.zeros(1)-1 #This will become bigger array. Set to -1 as check on whether it has been appended to yet
    vel_cluster_sizes = [] #number of nodes in clusters
    vel_cluster_areas = [] #number of surface-exposed nodes in clusters
    all_num_vel_clusters = np.zeros(1)-1 #This will become bigger array. Set to -1 as check on whether it has been appended to yet
    
    for subfolder in subfolder_list:
        print(subfolder)
        #Load data from reference lattice
        eq_traj = h5py.File(subfolder+'/equil/traj.h5')
        N, Lx, Ly, Lz, ref_lat = load_ref(eq_traj)
        adj_array = np.array(eq_traj['/connectivity/all_bonds'])

        #Extract neighbor lists for each node
        #nbs = np.zeros((N,12), dtype=int) #assumes fixed connectivity of 12
        nbs = []
        for i in range(N):
            nbs_i1 = np.unique(adj_array[adj_array[:,0]==i][:,1])
            nbs_i2 = np.unique(adj_array[adj_array[:,1]==i][:,0])
            nbs_i = np.concatenate([nbs_i1,nbs_i2])
            #check that this is 12 for the present fcc lattice
            if nbs_i.shape[0]!=12:
                print('Warning: node does not have normal fcc connectivity.')
            nbs.append(nbs_i)
        if len(nbs)!=N:
            print('Warning: neighbor list length does not equal # of nodes.')

        #Load trajectory
        traj = h5py.File(subfolder+'/prod/traj.h5')

        #Get magnitude of displacements and velocities
        disp_all = la.norm(np.array(traj['/particles/all/position/value'])-ref_lat,axis=2)
        vel_all = la.norm(np.array(traj['/particles/all/velocity/value']),axis=2)

        #Threshold based on input percentile
        disp_threshed = np.where(disp_all>disp_thresh,1,0)
        vel_threshed = np.where(vel_all>vel_thresh,1,0)

        #Identify "clusters" as contiguous collections of high displacements/velocities
        traj_length = disp_threshed.shape[0]
        num_disp_clusters = np.zeros(traj_length)
        num_vel_clusters = np.zeros(traj_length)

        for t in range(traj_length):
            if t%100==0:
                print('frame ', t)

            disp_cluster_id = get_clusters(disp_threshed[t,:], nbs)
            num_disp_clusters[t] = np.max(disp_cluster_id)
            disp_unique, disp_counts = np.unique(disp_cluster_id, return_counts=True)

            vel_cluster_id = get_clusters(vel_threshed[t,:], nbs)
            num_vel_clusters[t] = np.max(vel_cluster_id)
            vel_unique, vel_counts = np.unique(vel_cluster_id, return_counts=True)

            #Record volume (# of nodes in cluster)
            #and (TODO) surface area (# of exposed nodes in cluster)
            disp_cluster_sizes.append(disp_counts[1:]) #exclude the "0" cluster
            vel_cluster_sizes.append(vel_counts[1:]) #exclude the "0" cluster

        if np.any(all_num_disp_clusters == -1):
            all_num_disp_clusters = num_disp_clusters
        else:
            all_num_disp_clusters = np.vstack([all_num_disp_clusters, num_disp_clusters])

        if np.any(all_num_vel_clusters == -1):
            all_num_vel_clusters = num_vel_clusters
        else:
            all_num_vel_clusters = np.vstack([all_num_vel_clusters, num_vel_clusters])

    #Get histograms
    disp_cluster_sizes = np.concatenate(disp_cluster_sizes)
    all_num_disp_clusters = all_num_disp_clusters.flatten()

    disp_size_hist, disp_size_bin_edges = np.histogram(disp_cluster_sizes, np.arange(0,N+2,1)-0.5, density=True)
    disp_size_bins = (disp_size_bin_edges[:-1]+disp_size_bin_edges[1:])/2

    disp_num_hist, disp_num_bin_edges = np.histogram(all_num_disp_clusters, np.arange(0,N+2,1)-0.5, density=True)

    vel_cluster_sizes = np.concatenate(vel_cluster_sizes)
    all_num_vel_clusters = all_num_vel_clusters.flatten()

    vel_size_hist, vel_size_bin_edges = np.histogram(vel_cluster_sizes, np.arange(0,N+2,1)-0.5, density=True)
    vel_size_bins = (vel_size_bin_edges[:-1]+vel_size_bin_edges[1:])/2

    vel_num_hist, vel_num_bin_edges = np.histogram(all_num_vel_clusters, np.arange(0,N+2,1)-0.5, density=True)

    #Write data
    np.savetxt(in_folder + '/cluster_hist_disp_thresh=%f.txt' % percentile_thresh, np.c_[disp_size_bins,disp_size_hist,disp_num_hist], header='bin size num')
    np.savetxt(in_folder + '/cluster_hist_vel_thresh=%f.txt' % percentile_thresh, np.c_[vel_size_bins,vel_size_hist,vel_num_hist], header='bin size num') 

def load_ref(traj):

    N = traj['/particles/all/position/value'].shape[1]
    Lx = traj['/particles/all/box/edges'][0]
    Ly = traj['/particles/all/box/edges'][1]
    Lz = traj['/particles/all/box/edges'][2]
    ref_lat = traj['/particles/all/position/value'][0,:]

    return N, Lx, Ly, Lz, ref_lat

def extract_frame(traj, t, Lx, Ly, Lz):
    
    pos = traj['/particles/all/position/value'][t,:]
    vel = traj['/particles/all/velocity/value'][t,:]
    tstep = traj['/particles/all/step'][t]
    time = traj['/particles/all/time'][t]
    
    return tstep, time, pos, vel

@numba.jit(nopython=True)
def get_clusters(field, nbs):

    #Returns a numpy array specifying the index of the cluster
    #to which each node belongs

    N = field.shape[0]

    #print('Getting clusters...')
    cluster_id = np.zeros((N,),dtype=numba.int32)

    clusternumber = 0
    for i in range(N):
        if field[i]==1 and cluster_id[i]==0:
            clusternumber += 1
            cluster_id[i] = clusternumber
            harvest_cluster(clusternumber, i, cluster_id, field, nbs)

    return cluster_id

@numba.jit(nopython=True)
def harvest_cluster(clusternumber, ipart, cluster_id, field, nbs):

    for i in range(nbs[ipart].shape[0]):
        j = nbs[ipart][i]
        if field[j]==1 and cluster_id[j]==0:
            cluster_id[j] = clusternumber
            harvest_cluster(clusternumber, j, cluster_id, field, nbs)
        
def run_tests(in_folder):

    #Load data from reference lattice
    eq_traj = h5py.File(in_folder+'seed=1/equil/traj.h5')
    N, Lx, Ly, Lz, ref_lat = load_ref(eq_traj)
    adj_array = np.array(eq_traj['/connectivity/all_bonds'])

    #Extract neighbor lists for each node
    #nbs = np.zeros((N,12), dtype=int) #assumes fixed connectivity of 12
    nbs = []
    for i in range(N):
        nbs_i1 = np.unique(adj_array[adj_array[:,0]==i][:,1])
        nbs_i2 = np.unique(adj_array[adj_array[:,1]==i][:,0])
        nbs_i = np.concatenate([nbs_i1,nbs_i2])
        #check that this is 12 for the present fcc lattice
        if nbs_i.shape[0]!=12:
            print('Warning: node does not have normal fcc connectivity.')
        nbs.append(nbs_i)
    if len(nbs)!=N:
        print('Warning: neighbor list length does not equal # of nodes.')

    all_zeros = np.zeros(N)
    all_ones = np.ones(N)

    cluster_id_zeros = get_clusters(all_zeros, nbs)
    cluster_id_ones = get_clusters(all_ones, nbs)

    if(np.max(cluster_id_zeros)!=0):
        print('Test failed! Configuration with no nodes above threshold has cluster.')
        exit()
    if(np.max(cluster_id_ones)!=1):
        print('Test failed! Configuration with all nodes above threshold has no or more than one cluster.')
        exit()
    
    unique, counts = np.unique(cluster_id_ones, return_counts=True)
    if(counts[0]!=N):
        print('Test failed! Configuration with all nodes does not have cluster with N nodes.')
        exit()

    print('Tests passed!')
    
main()