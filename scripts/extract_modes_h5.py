#Get fourier components of displacement and velocity fields from h5md file

import numba
import sys
import os
import glob
import numpy as np
from math import sqrt
import numpy.linalg as la
import h5py

def main():

    h5_file = sys.argv[1] #trajectory file in h5md format
    ref_file = sys.argv[2] #trajectory file containing reference lattice
    out_folder = sys.argv[3] #where to put normal modes

    #Load data from reference lattice
    eq_traj = h5py.File(ref_file)
    N, Lx, Ly, Lz, a, ref_lat = load_ref(eq_traj)
    print(ref_lat)
    eq_traj.close()
    
    Nx = int(Lx/a)
    Ny = int(Ly/a)
    Nz = int(Lz/a)

    #reciprocal vectors
    q = get_q(N, Nx, Ny, Nz, a)
    indices = (np.rint(q/(2*np.pi)*Lx)).astype(int) #TODO: generalize to non-cubic box
    for i in range(q.shape[0]):
        print(indices[i,:])

    #Load trajectory
    traj = h5py.File(h5_file)
    pos_all = np.array(traj['/particles/all/position/value'])
    pos_all_unwrapped = get_pos_unwrapped(pos_all, np.array(traj['/particles/all/image/value']), Lx, Ly, Lz)
    vel_all = np.array(traj['/particles/all/velocity/value'])
    tsteps = np.array(traj['/particles/all/position/step'])
    times = np.array(traj['/particles/all/position/time'])
    traj.close()

    #Create output files
    uq_file = h5py.File(out_folder + '/uq.h5', 'w')
    vq_file = h5py.File(out_folder + '/vq.h5', 'w')

    uq_file.create_dataset('/parameters/a', data=a)
    uq_file.create_dataset('/parameters/edges', data=np.array([Lx,Ly,Lz]))
    uq_file.create_dataset('/indices', data=indices)

    vq_file.create_dataset('/parameters/a', data=a)
    vq_file.create_dataset('/parameters/edges', data=np.array([Lx,Ly,Lz]))
    vq_file.create_dataset('/indices', data=indices)

    #Loop through configuration in trajectory
    for t in range(pos_all.shape[0]):

        print('frame: ', t)
        
        #Load configuration
        tstep = tsteps[t]
        time = times[t]
        pos = pos_all[t,:,:]
        vel_r = vel_all[t,:,:]
        pos_unwrapped = pos_all_unwrapped[t,:,:]

        #tstep, time, pos, vel_r = extract_frame(traj, t, Lx, Ly, Lz)
            
        #compute displacement field
        disp_r = pos_unwrapped - ref_lat
        #Make sure to apply minimum image convention
        #disp_r = apply_min_image(disp_r, Lx, Ly, Lz)
        
        #Take discrete Fourier transform
        disp_q, vel_q = get_ft(disp_r, vel_r, ref_lat+np.array([Lx/2, Ly/2, Lz/2]), q, N)
        print('Done. Saving files...')

        #Now save the modes
        if t==0:
            uq_file.create_dataset('/data/step', data=np.array([tstep]), chunks=True, maxshape=(None,))
            uq_file.create_dataset('/data/time', data=np.array([time]), chunks=True, maxshape=(None,))
            uq_file.create_dataset('/data/value', data=np.array([disp_q]), chunks=True, maxshape=(None,disp_q.shape[0],disp_q.shape[1]))

            vq_file.create_dataset('/data/step', data=np.array([tstep]), chunks=True, maxshape=(None,))
            vq_file.create_dataset('/data/time', data=np.array([time]), chunks=True, maxshape=(None,))
            vq_file.create_dataset('/data/value', data=np.array([vel_q]), chunks=True, maxshape=(None,disp_q.shape[0],disp_q.shape[1]))
        else:
            print(uq_file['/data/step'].shape[0])
            uq_file['/data/step'].resize(uq_file['/data/step'].shape[0] + 1, axis=0)
            uq_file['/data/step'][-1] = tstep

            uq_file['/data/time'].resize(uq_file['/data/time'].shape[0] + 1, axis=0)
            uq_file['/data/time'][-1] = time

            uq_file['/data/value'].resize((uq_file['/data/value'].shape[0] + 1), axis=0)
            uq_file['/data/value'][-1] = disp_q

            vq_file['/data/step'].resize(vq_file['/data/step'].shape[0] + 1, axis=0)
            vq_file['/data/step'][-1] = tstep

            vq_file['/data/time'].resize(vq_file['/data/time'].shape[0] + 1, axis=0)
            vq_file['/data/time'][-1] = time

            vq_file['/data/value'].resize((vq_file['/data/value'].shape[0] + 1), axis=0)
            vq_file['/data/value'][-1] = vel_q

    #Display data
    #show_filtered_modes(q,disp_q)

    #test: try undoing FT
    #test_undo_ft(disp_r, disp_q, ref_lat, q, N)
    

def get_q(N, Nx, Ny, Nz, a):
    #Get allowed wavevectors (see Boykin et al. 2015)
    
    q = np.zeros((N,3))
    
    #define reciprocal lattice vectors
    b1 = (2*np.pi/a)*np.array((1.0, 1.0, -1.0))
    b2 = (2*np.pi/a)*np.array((1.0, -1.0, 1.0))
    b3 = (2*np.pi/a)*np.array((-1.0, 1.0, 1.0))
    b1_alt = (4*np.pi/a)*np.array((1.0, 0.0, 0.0))
    b2_alt = (4*np.pi/a)*np.array((0.0, 1.0, 0.0))
    b3_alt = (4*np.pi/a)*np.array((0.0, 0.0, 1.0))
    
    cnt = 0
    if Nx%2==1:
        nx_range = [-(Nx-1)//2 + i for i in range(Nx)]
    else:
        nx_range = [-(Nx-2)//2 + i for i in range(Nx)]
    if Ny%2==1:
        ny_range = [-(Ny-1)//2 + i for i in range(Ny)]
    else:
        ny_range = [-(Ny-2)//2 + i for i in range(Ny)]
    if Nz%2==1:
        nz_range = [-(Nz-1)//2 + i for i in range(Nz)]
    else:
        nz_range = [-(Nz-2)//2 + i for i in range(Nz)]
    for nx in nx_range:
        for ny in ny_range:
            for nz in nz_range:
                q[cnt,:] = (2*np.pi/a)*((nx/Nx)*np.array([1.0,0.0,0.0]) + (ny/Ny)*np.array([0.0,1.0,0.0]) + (nz/Nz)*np.array([0.0,0.0,1.0]))
                q[cnt+1,:] = q[cnt,:] + (2*np.pi)/a*np.array([1.0,0.0,0.0])
                q[cnt+2,:] = q[cnt,:] + (2*np.pi)/a*np.array([0.0,1.0,0.0])
                q[cnt+3,:] = q[cnt,:] + (2*np.pi)/a*np.array([0.0,0.0,1.0])
                cnt += 4
                
    #Shift vectors outside the first BZ back inside
    print(q)
    print('shifting wavevectors to first BZ...')
    for i in range(N):   
        if not in_first_bz(q[i], a):
            #print('vector %d is outside of 1st Brillouin zone' % i)
            reciprocal_list = [-b1_alt, b1_alt, -b2_alt, b2_alt, -b3_alt, b3_alt]
            q[i] = shift_to_first_bz(q[i],a,reciprocal_list)
    #Test that vectors are actually in 1st BZ
    for i in range(N):
        if not in_first_bz(q[i], a):
            print('Error: vector %d is outside of 1st Brillouin zone' % i)
            exit()
    #Check that wavevectors are unique
    if np.unique(q,axis=0, return_counts=True)[0].shape[0]!=N:
        print('Error: not all wavevectors are unique!')
        
    print('Got allowed q.')
    return q

def shift_to_first_bz(q, a, mylist, multiplier=1):
    for b in mylist:
        b = multiplier*b
        q_trial = q + b
        if in_first_bz(q_trial, a):
            #', q, 'to', q_trial)
            q = q_trial
            break
        if (b==mylist[-1]).all():
            #print('Going to multiple %d of reciprocal vectors...' % (multiplier+1))
            return shift_to_first_bz(q, a, mylist, multiplier+1)
    return q

def in_first_bz(q, a): 
    #args: vector to be tested, lattice constant
    
    X_pts = [np.array([-2*np.pi/a, 0, 0]),
             np.array([2*np.pi/a, 0, 0]),
             np.array([0, -2*np.pi/a, 0]),
             np.array([0, 2*np.pi/a, 0]),
             np.array([0, 0, -2*np.pi/a]),
             np.array([0, 0, 2*np.pi/a])]
    L_pts = [np.array([-2*np.pi/a,-2*np.pi/a,-2*np.pi/a]),
             np.array([2*np.pi/a,-2*np.pi/a,-2*np.pi/a]),
             np.array([-2*np.pi/a,2*np.pi/a,-2*np.pi/a]),
             np.array([-2*np.pi/a,-2*np.pi/a,2*np.pi/a]),
             np.array([2*np.pi/a,2*np.pi/a,-2*np.pi/a]),
             np.array([2*np.pi/a,-2*np.pi/a,2*np.pi/a]),
             np.array([-2*np.pi/a,2*np.pi/a,2*np.pi/a]),
             np.array([2*np.pi/a,2*np.pi/a,2*np.pi/a])]
    
    for q0 in X_pts:
        n = q0/la.norm(q0)
        z = np.dot(q-q0,n)
        if z>0:
            return False
        
    for q0 in L_pts:
        n = q0/la.norm(q0)
        z = np.dot(q-q0,n)
        if z>0:
            return False
    
    return True
    
#TEST: filter out "noise" to check for "signal"
def show_filtered_modes(q,uq):

    test_q = q[np.abs(uq[:,0])>1e-12]
    test_uq = uq[np.abs(uq[:,0])>1e-12]
    for i in range(test_q.shape[0]):
        print(i, test_q[i,:], test_uq[i,:])

#TEST: make sure undoing FT gives expected result
def test_undo_ft(disp_r, disp_q, ref_lat, q, N):

    disp_r_test = np.zeros((N,3), dtype=complex)
    for i in range(N):
        for j in range(N):
            disp_r_test[i,:] += (1.0/N)*disp_q[j,:]*np.exp(-1j*np.dot(ref_lat[i,:], q[j,:]))
    diff = disp_r_test-disp_r.astype(complex)
    if (np.abs(diff)>1e-13).any():
        print('Error: fourier transforms done improperly.')
        exit()

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
    #print(rdiff)
    rdiff = np.where(rdiff>arr1, rdiff-arr3, rdiff)
    rdiff = np.where(rdiff<arr2, rdiff+arr3, rdiff)
    return rdiff

def load_ref(traj):

    N = int(traj['/particles/all/position/value'].shape[1])
    Lx = float(traj['/particles/all/box/edges'][0])
    Ly = float(traj['/particles/all/box/edges'][1])
    Lz = float(traj['/particles/all/box/edges'][2])
    a = traj['/particles/all/box/lattice_constant']
    a = np.array(a)
    print(a.dtype)
    ref_lat = np.array(traj['/particles/all/position/value'][0,:])

    return N, Lx, Ly, Lz, a, ref_lat

def extract_frame(traj, t, Lx, Ly, Lz):
    
    pos = np.array(traj['/particles/all/position/value'][t,:])
    vel = np.array(traj['/particles/all/velocity/value'][t,:])
    tstep = int(traj['/particles/all/position/step'][t])
    time = float(traj['/particles/all/position/time'][t])
    
    return tstep, time, pos, vel
    
@numba.jit(nopython=True)
def get_ft(disp_r, vel_r,ref_lat, q, N):
    print('Getting FT...')
    disp_q = np.zeros((N,3), dtype=np.complex128)
    vel_q = np.zeros((N,3), dtype=np.complex128)
    for i in range(N):
        for j in range(N):
            disp_q[i,:] += disp_r[j,:]*np.exp(1j*np.dot(ref_lat[j,:], q[i,:]))
            vel_q[i,:] += vel_r[j,:]*np.exp(1j*np.dot(ref_lat[j,:], q[i,:]))
    return disp_q, vel_q
        
@numba.jit(nopython=True)
def get_pos_unwrapped(pos, image, Lx, Ly, Lz):

    shift = np.zeros(image.shape)
    shift[:,:,0] = image[:,:,0]*Lx
    shift[:,:,1] = image[:,:,1]*Ly
    shift[:,:,2] = image[:,:,2]*Lz

    pos_unwrapped = pos + shift
    
    return pos_unwrapped

main()
