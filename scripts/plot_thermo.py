#Plot thermodynamic quantities (energy, ...)

import numpy as np
import pylab as plt
import sys

#TODO: make these more specific
dat_file = sys.argv[1] #should be a .dat file
out_dir = sys.argv[2] #should be a directory

data = np.loadtxt(dat_file, skiprows=1) #ignore header

#TODO: read these in from file
N = 500
kT = 0.2
K = 20.0

Eapprox = (3.0/2.0)*N*kT

#Energy vs. time
#TODO: read in parameters from file
fig = plt.figure()
plt.plot(data[:,0],data[:,1],label=r'Langevin, $\gamma=10$, $dt=1 \times 10^{-3}$')
plt.xlabel('timestep')
plt.ylabel(r'$E$')
plt.title(r'$k_BT/\epsilon=%.02f,$ $N=%d,$ $Kl_0^2/\epsilon=%.01f$' % (kT, N, K))
#plt.xlim([0,10])
#plt.ylim([-100,100])
plt.axhline(y=Eapprox,color='black',label='Equipartition')
plt.legend(fontsize=10)
plt.savefig(out_dir + '/' + 'energy_vs_t.pdf')
plt.show()