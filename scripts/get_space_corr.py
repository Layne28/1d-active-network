#Compute spatial correlation functions

import numba
import sys
import os
import glob
import numpy as np
from math import sqrt
import numpy.linalg as la

def main():

    in_folder = sys.argv[1] #folder containing dump files
    ref_file = sys.argv[2] #dump file for reference lattice
    out_folder = sys.argv[3]
    #thetype = sys.argv[4] #displacement or velocity

    '''
    if thetype=='disp':
        name='u'
        print('Getting time correlation of displacement modes.')
    elif thetype=='vel':
        name='v'
        print('Getting time correlation of velocity modes.')
    else:
        print('Error: type of mode not recognized.')
        exit()
    '''

    #Load data from reference lattice
    time, N, Lx, Ly, Lz, map_atom_to_pos = load_ref_dump(ref_file)
    ref_lat = np.array(list(map_atom_to_pos.values()))
    
    #Map positions to lattice indices
    Nx = round((N/4)**(1.0/3.0)) #Assumes cube! TODO: generalize to arbitrary box shape
    Ny = Nx
    Nz = Nx
    print(Nx)
    a = Lx/Nx #lattice constant

    #test get_min_disp
    test1 = np.array([1.1,0.3,-1.2])
    test2 = np.array([0,0,0])
    print('testing get_min_disp...')
    print(get_min_disp(test1, test2, 2,2,2))

    cuu_r = {}
    cvv_r = {}
    count_r = {}
    #Get possible separations
    for i in range(N):
        for j in range(i,N):
            sep = get_min_disp(ref_lat[j,:],ref_lat[i,:],Lx,Ly,Lz)
            cuu_r[(format(sep[0], '.6f'),format(sep[1], '.6f'),format(sep[2], '.6f'))] = np.zeros((3,3))
            cvv_r[(format(sep[0], '.6f'),format(sep[1], '.6f'),format(sep[2], '.6f'))] = np.zeros((3,3))
            count_r[(format(sep[0], '.6f'),format(sep[1], '.6f'),format(sep[2], '.6f'))] = 0
    for i in range(N):
        for j in range(i,N):
            sep = get_min_disp(ref_lat[j,:],ref_lat[i,:],Lx,Ly,Lz)
            count_r[(format(sep[0], '.6f'),format(sep[1], '.6f'),format(sep[2], '.6f'))] += 1

    if(os.path.exists(in_folder)):
        mylist = []
        for thefile in sorted(sorted(glob.glob(in_folder + "/*.dump")), key=len):
            #print(thefile)
            mylist.append(thefile)
        print('Loaded files.')
    else:
        print('Input folder does not exist. Exiting.')
        exit()

    if not os.path.exists(out_folder):
        print('Output folder does not exist. Creating it now.')
        os.makedirs(out_folder)

    cuu = np.zeros((N*(N+1)//2,3,3))
    cvv = np.zeros((N*(N+1)//2,3,3))
    #Loop through files
    for thefile in mylist:

        print(thefile)
        
        #Load configuration
        tstep, pos, vel_r = load_dump(thefile)
            
        #compute displacement field
        disp_r = pos - ref_lat #Need to apply minimum image convention
        disp_r = apply_min_image(disp_r, Lx, Ly, Lz)

        sep_arr, cuu_arr, cvv_arr = get_corr(disp_r, vel_r, ref_lat, N, Lx, Ly, Lz)

        cuu += cuu_arr
        cvv += cvv_arr


        '''
        for i in range(N):
            for j in range(i,N):
                sep = get_min_disp(ref_lat[j,:],ref_lat[i,:],Lx,Ly,Lz)
                cuu_r[(format(sep[0], '.6f'),format(sep[1], '.6f'),format(sep[2], '.6f'))] += disp_r[i,:,None]*disp_r[j,None,:]
                cvv_r[(format(sep[0], '.6f'),format(sep[1], '.6f'),format(sep[2], '.6f'))] += vel_r[i,:,None]*vel_r[j,None,:]
        '''
    for i in range(sep_arr.shape[0]):
        cuu_r[(format(sep_arr[i,0], '.6f'),format(sep_arr[i,1], '.6f'),format(sep_arr[i,2], '.6f'))] += cuu[i,:,:]
        cvv_r[(format(sep_arr[i,0], '.6f'),format(sep_arr[i,1], '.6f'),format(sep_arr[i,2], '.6f'))] += cvv[i,:,:]
    print(len(cvv_r.keys()))
    print(len(cuu_r.keys()))
    for key in count_r.keys():
        cuu_r[key] /= (len(mylist)*count_r[key])
        cvv_r[key] /= (len(mylist)*count_r[key])
    print(np.array(count_r.keys()))

    pos_arr = np.zeros((len(count_r.keys()),4))
    print(pos_arr.shape)
    ucorr_arr = np.zeros((len(count_r.keys()),3,3))
    vcorr_arr = np.zeros((len(count_r.keys()),3,3))

    cnt = 0
    for key in count_r.keys():
        pos_arr[cnt,0] = float(key[0])
        pos_arr[cnt,1] = float(key[1])
        pos_arr[cnt,2] = float(key[2])
        pos_arr[cnt,3] = np.sqrt(pos_arr[cnt,0]**2+pos_arr[cnt,1]**2+pos_arr[cnt,2]**2)
        ucorr_arr[cnt,:,:] = cuu_r[key]
        vcorr_arr[cnt,:,:] = cvv_r[key]
        #print(pos_arr[cnt,:])
        cnt += 1
    newinds = pos_arr[:,-1].argsort()
    pos_arr = pos_arr[newinds]
    ucorr_arr = ucorr_arr[newinds]
    vcorr_arr = vcorr_arr[newinds]
    print(ucorr_arr)
    print(vcorr_arr)

    np.savez(out_folder + '/u_r_corr.npz', pos=pos_arr,corr=ucorr_arr)
    np.savez(out_folder + '/v_r_corr.npz', pos=pos_arr,corr=vcorr_arr)

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
def get_corr_alt(cuu_r, cvv_r, disp_r, vel_r, ref_lat, N, Lx, Ly, Lz):
    for i in range(N):
        for j in range(i,N):
            sep = get_min_disp(ref_lat[j,:],ref_lat[i,:],Lx,Ly,Lz)
            cuu_r[(format(sep[0], '.6f'),format(sep[1], '.6f'),format(sep[2], '.6f'))] += disp_r[i,:,None]*disp_r[j,None,:]
            cvv_r[(format(sep[0], '.6f'),format(sep[1], '.6f'),format(sep[2], '.6f'))] += vel_r[i,:,None]*vel_r[j,None,:]

@numba.jit(nopython=True) 
def get_corr(disp_r, vel_r, ref_lat, N, Lx, Ly, Lz):
    sep_arr = np.zeros((N*(N+1)//2,3))
    cuu_arr = np.zeros((N*(N+1)//2,3,3))
    cvv_arr = np.zeros((N*(N+1)//2,3,3))
    
    cnt=0
    for i in range(N):
        for j in range(i,N):
            sep_arr[cnt,:] = get_min_disp(ref_lat[j,:],ref_lat[i,:],Lx,Ly,Lz)
            for n1 in range(3):
                for n2 in range(3):
                    cuu_arr[cnt,n1,n2] += disp_r[i,n1]*disp_r[j,n2]
                    cvv_arr[cnt,n1,n2] += vel_r[i,n1]*vel_r[j,n2]
            cnt += 1

    return sep_arr, cuu_arr, cvv_arr

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
        
main()
