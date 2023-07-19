#Get a mode trajectory for each wavevector

import sys
import os
import glob
import numpy as np
from math import sqrt
import numpy.linalg as la

def main():

    in_folder = sys.argv[1]
    out_folder = sys.argv[2]

    if(os.path.exists(in_folder)):
        mylist = []
        for thefile in sorted(sorted(glob.glob(in_folder + "/uq*.txt")), key=len):
            print(thefile)
            mylist.append(thefile)
    else:
        print('Input folder does not exist. Exiting.')
        exit()

    if not os.path.exists(out_folder):
        print('Output folder does not exist. Creating it now.')
        os.makedirs(out_folder)

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
    print('Loaded data.')

    #Get indices of allowed qs
    #assumes lattice constant a=sqrt(2) and Nx=5
    qs = uqs[0][:,:3]
    indices = np.rint(np.real(qs)*(sqrt(2.0)/(2*np.pi))*5).astype(int)
    print(indices)

    for i in range(indices.shape[0]):
        n1 = indices[i,0]
        n2 = indices[i,1]
        n3 = indices[i,2]
        uqmat = np.array([uq[i,3:] for uq in uqs])
        np.savez(out_folder + '/uq_traj_nx=%d_ny=%d_nz=%d.npz' % (n1, n2, n3), time=times, traj=uqmat)

main()