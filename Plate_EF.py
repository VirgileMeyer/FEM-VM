import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from matplotlib.mlab import griddata
from mpl_toolkits.mplot3d import Axes3D
from plate_value import *
import timeit
import pickle


np.set_printoptions(threshold=np.nan)

omega=2*np.pi*freq

#element
a=0
b=0
element=np.zeros(((nx-1)*(ny-1),13))

for i in range(0,nx-2+1):
    for j in range(0,ny-2+1):
        a=3*i+3*j*nx
        p11=(1+a)+b
        p12=(1+a)+b+1
        p13=(1+a)+b+2
        p21=(3*nx+1+a)+b
        p22=(3*nx+1+a)+b+1
        p23=(3*nx+1+a)+b+2
        p31=(3*nx+2+a)+b+2
        p32=(3*nx+2+a)+b+3
        p33=(3*nx+2+a)+b+4
        p41=(2+a)+b+2
        p42=(2+a)+b+3
        p43=(2+a)+b+4
        elem=1+i+j*(nx-1)
        element[i+j*(nx-1)]=[elem,p11,p12,p13,p41,p42,p43,p31,p32,p33,p21,p22,p23]
        
#assembling

kp=np.full((3*nx*ny,3*nx*ny),0.)
mp=np.full((3*nx*ny,3*nx*ny),0.)

for el in range (0,(nx-1)*(ny-1)):
    for a in range (0,12,3):
         for b in range (0,12,3):
             kp[int(element[el,b+1]-1),int(element[el,a+1]-1)]=kp[int(element[el,b+1]-1),int(element[el,a+1]-1)]+kep[b,a]
             kp[int(element[el,b+2]-1),int(element[el,a+2]-1)]=kp[int(element[el,b+2]-1),int(element[el,a+2]-1)]+kep[b+1,a+1]
             kp[int(element[el,b+3]-1),int(element[el,a+3]-1)]=kp[int(element[el,b+3]-1),int(element[el,a+3]-1)]+kep[b+2,a+2]
             kp[int(element[el,b+1]-1),int(element[el,a+2]-1)]=kp[int(element[el,b+1]-1),int(element[el,a+2]-1)]+kep[b,a+1]
             kp[int(element[el,b+2]-1),int(element[el,a+1]-1)]=kp[int(element[el,b+2]-1),int(element[el,a+1]-1)]+kep[b+1,a]
             kp[int(element[el,b+2]-1),int(element[el,a+3]-1)]=kp[int(element[el,b+2]-1),int(element[el,a+3]-1)]+kep[b+1,a+2]
             kp[int(element[el,b+3]-1),int(element[el,a+2]-1)]=kp[int(element[el,b+3]-1),int(element[el,a+2]-1)]+kep[b+2,a+1]
             kp[int(element[el,b+3]-1),int(element[el,a+1]-1)]=kp[int(element[el,b+3]-1),int(element[el,a+1]-1)]+kep[b+2,a]
             kp[int(element[el,b+1]-1),int(element[el,a+3]-1)]=kp[int(element[el,b+1]-1),int(element[el,a+3]-1)]+kep[b,a+2]

for el in range (0,(nx-1)*(ny-1)):
    for a in range (0,12,3):
         for b in range (0,12,3):
             mp[int(element[el,b+1]-1),int(element[el,a+1]-1)]=mp[int(element[el,b+1]-1),int(element[el,a+1]-1)]+mep[b,a]
             mp[int(element[el,b+2]-1),int(element[el,a+2]-1)]=mp[int(element[el,b+2]-1),int(element[el,a+2]-1)]+mep[b+1,a+1]
             mp[int(element[el,b+3]-1),int(element[el,a+3]-1)]=mp[int(element[el,b+3]-1),int(element[el,a+3]-1)]+mep[b+2,a+2]
             mp[int(element[el,b+1]-1),int(element[el,a+2]-1)]=mp[int(element[el,b+1]-1),int(element[el,a+2]-1)]+mep[b,a+1]
             mp[int(element[el,b+2]-1),int(element[el,a+1]-1)]=mp[int(element[el,b+2]-1),int(element[el,a+1]-1)]+mep[b+1,a]
             mp[int(element[el,b+2]-1),int(element[el,a+3]-1)]=mp[int(element[el,b+2]-1),int(element[el,a+3]-1)]+mep[b+1,a+2]
             mp[int(element[el,b+3]-1),int(element[el,a+2]-1)]=mp[int(element[el,b+3]-1),int(element[el,a+2]-1)]+mep[b+2,a+1]
             mp[int(element[el,b+3]-1),int(element[el,a+1]-1)]=mp[int(element[el,b+3]-1),int(element[el,a+1]-1)]+mep[b+2,a]
             mp[int(element[el,b+1]-1),int(element[el,a+3]-1)]=mp[int(element[el,b+1]-1),int(element[el,a+3]-1)]+mep[b,a+2]
             

#global matrix
globplate=Di*kp-rost*omega**2*h*mp

# boundary conditions CCCC
globplate1=globplate
for i in range (0,3*nx-1 +1):
    globplate1=np.delete(globplate1,3*nx*ny-i-1,0)
    globplate1=np.delete(globplate1,3*nx*ny-i-1,1)

for i in range (0,3*nx*(ny-3) +1,3*nx):
    globplate1=np.delete(globplate1,3*(ny-1)*nx-i-1,0)
    globplate1=np.delete(globplate1,3*(ny-1)*nx-i-1,1)
    
    globplate1=np.delete(globplate1,3*(ny-1)*nx-i-1 -1,0)
    globplate1=np.delete(globplate1,3*(ny-1)*nx-i-1 -1,1)

    globplate1=np.delete(globplate1,3*(ny-1)*nx-i-2 -1,0)
    globplate1=np.delete(globplate1,3*(ny-1)*nx-i-2 -1,1)

    globplate1=np.delete(globplate1,3*(ny-2)*nx+1-i+2 -1,0)
    globplate1=np.delete(globplate1,3*(ny-2)*nx+1-i+2 -1,1)

    globplate1=np.delete(globplate1,3*(ny-2)*nx+1-i+1 -1,0)
    globplate1=np.delete(globplate1,3*(ny-2)*nx+1-i+1 -1,1)

    globplate1=np.delete(globplate1,3*(ny-2)*nx+1-i -1,0)
    globplate1=np.delete(globplate1,3*(ny-2)*nx+1-i -1,1)

for i in range(0,3*nx-1 +1):
    globplate1=np.delete(globplate1,3*nx-i -1,0)
    globplate1=np.delete(globplate1,3*nx-i -1,1)



# "inverted" = (Di*kp - rhost*h*omega*2*mp)**-1

platinv=np.linalg.inv(globplate1)

platinvdepl=np.full((nx*ny-2*nx-2*(ny-2),nx*ny-2*nx-2*(ny-2)),0.)

for i in range(0,nx*ny-2*nx-2*(ny-2)-1+1):
    for j in range(0,nx*ny-2*nx-2*(ny-2)-1+1):
        platinvdepl[i,j]=platinv[3*i,3*j]

# displacement is obtained by computing "inverted" * f*exp(i*omega*t)





# print(platinvdepl)

xs=0.05
ys=0.25
node=0

nodesource=None

# loop start at 1 to take into account CCCC
for i in range (1, nx-1):
    for j in range (1, ny-1):
        node=node+1
        x=dx*(i-1)
        y=dy*(j-1)
        if abs(x-xs)<dx:
            if abs(y-ys)<dy:
                nodesource=node
                
if nodesource == None:
    print("Source error")
else:
    print("xs=",xs," ","ys=",ys," ","nodesource=",nodesource)


# right hand side force
fb=np.zeros(nx*ny-2*nx-2*(ny-2))
fb[nodesource-1]=1
print("nodesource =",nodesource)

print("total nodes minus BC",nx*ny-2*nx-2*(ny-2))

deplfield=np.dot(platinvdepl,fb)

node=0
plotfile=np.zeros((nx*ny-2*nx-2*(ny-2),3))               
for i in range(0,nx-2):
    for j in range(0,ny-2):
        x=dx*i
        y=dy*j
        plotfile[node]=[x,y,deplfield[node]]
        node=node+1
# print(deplfield)

# To save data for further plot
pickle.dump(plotfile, open('deplfield.p', 'wb'))


x=plotfile[:,0]
y=plotfile[:,1]
z=plotfile[:,2]

# 3D plot
fig = plt.figure()
ax = Axes3D(fig)
surf = ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.1)
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()


# displacement field along y-axis
i0=5
deplfdaxis=np.full(ny-2,0.)
yaxis=np.full(ny-2,0.)
for i in range (0,ny-2):
    j=i*(nx-2)
    yaxis[i]=dy*i
    deplfdaxis[i]=deplfield[i0+j]
                    
