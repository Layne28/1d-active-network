import numba
import sys
import os
import glob
import numpy as np
from math import sqrt
import numpy.linalg as la
import h5py

def main():

    in_folder = sys.argv[1]

    print(in_folder)

    #Load trajectory
    uq = h5py.File(in_folder + '/uq.h5')
    vq = h5py.File(in_folder + '/vq.h5')

    print('uq size: ', uq['/data/value'].shape[0])
    print('vq size: ', vq['/data/value'].shape[0])

    with open('check_size.txt', 'a') as f:
        f.write(in_folder+'\n')
        f.write('uq size: ' + str(uq['/data/value'].shape[0]) + '\n')
        f.write('vq size: ' + str(vq['/data/value'].shape[0]) + '\n')

main()
