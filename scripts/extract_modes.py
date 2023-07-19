#Get fourier components of displacement and velocity fields from lammps dump file

import numba
import sys
import os
import glob
import numpy as np
from math import sqrt
import numpy.linalg as la

def main():

    #Input. If "pos_file" ends with ".dump" this will take a single FT. Otherwise it will interpret "pos_file" as a folder containing files and try looping through all the ".dump" files in that folder
    ref_file = sys.argv[1] #dump file for reference lattice
    pos_file = sys.argv[2] #either a single network configuration or a folder containing many
    out_folder = sys.argv[3] #where to put normal modes

    #Load data from reference lattice
    time, N, Lx, Ly, Lz, map_atom_to_pos = load_ref_dump(ref_file)
    ref_lat = np.array(list(map_atom_to_pos.values()))
    print(ref_lat)
    
    #Map positions to lattice indices
    Nx = round((N/4)**(1.0/3.0)) #Assumes cube! TODO: generalize to arbitrary box shape
    Ny = Nx
    Nz = Nx
    print(Nx)
    a = Lx/Nx #lattice constant

    #reciprocal vectors
    q = get_q(N, Nx, Ny, Nz, a)
    indices = (np.rint(q*sqrt(2.0)/(2*np.pi)*Nx)).astype(int)
    for i in range(q.shape[0]):
        print(indices[i,:])
    #exit()

    #check whether this is a single file or directory
    if pos_file.endswith('.dump'):
        mylist = [pos_file]
    elif(os.path.exists(pos_file)):
        mylist = []
        for thefile in sorted(sorted(glob.glob(pos_file + "/*.dump")), key=len):
            print(thefile)
            mylist.append(thefile)
    else:
        print('Error: argument is neither dump file nor directory. Exiting.')
        exit()

    #Loop through files
    for thefile in mylist:

        print(thefile)
        
        #Load configuration
        tstep, pos, vel_r = load_dump(thefile)
            
        #compute displacement field
        disp_r = pos - ref_lat
        #Make sure to apply minimum image convention
        disp_r = apply_min_image(disp_r, Lx, Ly, Lz)
        
        #Take discrete Fourier transform
        disp_q, vel_q = get_ft(disp_r, vel_r, ref_lat, q, N)
        print('Done. Saving files...')

        np.savetxt(out_folder + '/uq_%s.txt' % (thefile.split('/')[-1]).split('.')[0], np.c_[q,disp_q], header='timestep: %d' % tstep)
        np.savetxt(out_folder + '/vq_%s.txt' % (thefile.split('/')[-1]).split('.')[0], np.c_[q,vel_q], header='timestep: %d' % tstep)

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

def load_ref_dump(myfile):
    
    if not myfile.endswith('.dump'):
        print('Error: must provide dump file.')
        exit()
    
    with open(myfile) as f:
        lines = f.readlines()
        
    time = float(lines[1])
    N = int(lines[3])
    Lx = float((lines[5].split())[-1])-float((lines[5].split())[0])
    Ly = float((lines[6].split())[-1])-float((lines[6].split())[0])
    Lz = float((lines[7].split())[-1])-float((lines[7].split())[0])
    
    #Get map from atom index to lattice vector
    mymap = {}
    for line in lines[9:]:
        index = int((line.split())[0])
        x = float((line.split())[3])
        y = float((line.split())[4])
        z = float((line.split())[5])
        mymap[index] = np.array((x, y, z))+np.array([Lx/2, Ly/2, Lz/2])

    return time, N, Lx, Ly, Lz, mymap

def load_dump(myfile):
    
    if not myfile.endswith('.dump'):
        print('Error: must provide dump file.')
        exit()
    
    with open(myfile) as f:
        lines = f.readlines()
        
    #TODO: check that these agree with reference lattice
    tstep = int(lines[1])
    N = int(lines[3])
    Lx = float((lines[5].split())[-1])-float((lines[5].split())[0])
    Ly = float((lines[6].split())[-1])-float((lines[6].split())[0])
    Lz = float((lines[7].split())[-1])-float((lines[7].split())[0])
    pos = np.zeros((N,3))
    vel = np.zeros((N,3))
    
    for line in lines[9:]:
        index = int((line.split())[0])
        x = float((line.split())[3])
        y = float((line.split())[4])
        z = float((line.split())[5])
        vx = float((line.split())[6])
        vy = float((line.split())[7])
        vz = float((line.split())[8])
        pos[index,0] = x
        pos[index,1] = y
        pos[index,2] = z
        pos[index,:] += np.array([Lx/2, Ly/2, Lz/2])
        vel[index,0] = vx
        vel[index,1] = vy
        vel[index,2] = vz
    
    return tstep, pos, vel
    
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
        
main()
