#Create a dump file with a deformation
#at a particular wavevector for testing

import sys
import numpy as np
from math import sqrt
import numpy.linalg as la

def main():
    #Load data from reference lattice
    time, N, Lx, Ly, Lz, pos = load_ref_dump(sys.argv[1])
    thefolder = sys.argv[1][:(sys.argv[1].rindex('/')+1)]

    print(N)

    #Apply long-wavelength deformation along x axis
    for i in range(N):
        pos[i,0] += 0.2*np.cos(2*np.pi*pos[i,0]/Lx)

    #Save test configuration
    save_dump(thefolder + 'test_deformation.dump', time, N, Lx, Ly, Lz, pos)

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
    
    #Get positions
    pos = np.zeros((N,3))
    for line in lines[9:]:
        index = int((line.split())[0])
        x = float((line.split())[3])
        y = float((line.split())[4])
        z = float((line.split())[5])
        pos[index,:] = np.array((x, y, z))

    return time, N, Lx, Ly, Lz, pos

def load_dump(myfile):
    
    if not myfile.endswith('.dump'):
        print('Error: must provide dump file.')
        exit()
    
    with open(myfile) as f:
        lines = f.readlines()
        
    #TODO: check that these agree with reference lattice
    N = int(lines[3])
    Lx = float((lines[5].split())[-1])-float((lines[5].split())[0])
    Ly = float((lines[6].split())[-1])-float((lines[6].split())[0])
    Lz = float((lines[7].split())[-1])-float((lines[7].split())[0])
    pos = np.zeros((N,3))
    
    for line in lines[9:]:
        index = int((line.split())[0])
        x = float((line.split())[3])
        y = float((line.split())[4])
        z = float((line.split())[5])
        pos[index,0] = x
        pos[index,1] = y
        pos[index,2] = z
        pos[index,:] += np.array([Lx/2, Ly/2, Lz/2])
    
    return pos
    
def save_dump(myfile, t, N, Lx, Ly, Lz, pos):

    with open(myfile, 'w') as f:
        f.write('ITEM: TIMESTEP\n')
        f.write('%d\n' % t)
        f.write('ITEM: NUMBER OF ATOMS\n')
        f.write('%d\n' % N)
        f.write('ITEM: BOX BOUNDS pp pp pp\n')
        f.write('%.15e %.15e\n' % (-Lx/2, Lx/2))
        f.write('%.15e %.15e\n' % (-Ly/2, Ly/2))
        f.write('%.15e %.15e\n' % (-Lz/2, Lz/2))
        f.write('ITEM: ATOMS id mol type x y z\n')
        for i in range(N):
            f.write('%d 1 1 %.15e %.15e %.15e\n' % (i, pos[i,0], pos[i,1], pos[i,2]))

        
main()