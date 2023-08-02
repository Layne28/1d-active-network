#Get strain histogram from collection of h5 files

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import os
import glob



def main():

    in_folder = sys.argv[1]

    data = np.load(in_folder + '/strain_hist.npz')
    e = {name: data[name] for name in data}
    print(e)

    fig = plt.figure()
    nchunks = data['nchunks']
    colors = cm.viridis(np.linspace(0,1,nchunks))
    for n in range(nchunks):
        plt.plot(data['bins_%d' % n], data['hist_%d' % n], color=colors[n])
    plt.xlabel(r'$s$')
    plt.ylabel(r'$P(s)$')
    plt.savefig('strain_histo.png',dpi=300,bbox_inches='tight')
    plt.show()

main()