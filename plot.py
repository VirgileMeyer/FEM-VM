import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from matplotlib.mlab import griddata
from mpl_toolkits.mplot3d import Axes3D


plotfile=pickle.load(open('deplfield.p', 'rb'))


x=plotfile[:,0]
y=plotfile[:,1]
z=plotfile[:,2]


#fig = plt.figure()
#ax = Axes3D(fig)
#surf = ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.1)
#fig.colorbar(surf, shrink=0.5, aspect=5)



plt.figure()
plt.tricontour(x, y, z)

plt.show()




