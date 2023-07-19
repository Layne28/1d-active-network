#Plot thermodynamic quantities (energy, ...)

import numpy as np
import pylab as plt
import sys

in_file = sys.argv[1]

data = np.load(in_file)
print(data['corr'].shape)

dot = np.trace(data['corr'],axis1=1,axis2=2)
print(dot.shape)

fig = plt.figure()
plt.plot(data['time'],dot)
plt.show()
