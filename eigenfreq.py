import numpy as np
from Plate_EF import kp,mp
from plate_value import nx,ny,Di,rost,h


# compute the eigenfrequencies for a plate CCCC

kp1=kp
mp1=mp
for i in range (0,3*nx-1 +1):
    kp1=np.delete(kp1,3*nx*ny-i-1,0)
    kp1=np.delete(kp1,3*nx*ny-i-1,1)
for i in range (0,3*nx*(ny-3) +1,3*nx):
    kp1=np.delete(kp1,3*(ny-1)*nx-i-1,0)
    kp1=np.delete(kp1,3*(ny-1)*nx-i-1,1)
    kp1=np.delete(kp1,3*(ny-1)*nx-i-1 -1,0)
    kp1=np.delete(kp1,3*(ny-1)*nx-i-1 -1,1)
    kp1=np.delete(kp1,3*(ny-1)*nx-i-2 -1,0)
    kp1=np.delete(kp1,3*(ny-1)*nx-i-2 -1,1)
    kp1=np.delete(kp1,3*(ny-2)*nx+1-i+2 -1,0)
    kp1=np.delete(kp1,3*(ny-2)*nx+1-i+2 -1,1)
    kp1=np.delete(kp1,3*(ny-2)*nx+1-i+1 -1,0)
    kp1=np.delete(kp1,3*(ny-2)*nx+1-i+1 -1,1)
    kp1=np.delete(kp1,3*(ny-2)*nx+1-i -1,0)
    kp1=np.delete(kp1,3*(ny-2)*nx+1-i -1,1)
for i in range(0,3*nx-1 +1):
    kp1=np.delete(kp1,3*nx-i -1,0)
    kp1=np.delete(kp1,3*nx-i -1,1)

for i in range (0,3*nx-1 +1):
    mp1=np.delete(mp1,3*nx*ny-i-1,0)
    mp1=np.delete(mp1,3*nx*ny-i-1,1)
for i in range (0,3*nx*(ny-3) +1,3*nx):
    mp1=np.delete(mp1,3*(ny-1)*nx-i-1,0)
    mp1=np.delete(mp1,3*(ny-1)*nx-i-1,1)
    mp1=np.delete(mp1,3*(ny-1)*nx-i-1 -1,0)
    mp1=np.delete(mp1,3*(ny-1)*nx-i-1 -1,1)
    mp1=np.delete(mp1,3*(ny-1)*nx-i-2 -1,0)
    mp1=np.delete(mp1,3*(ny-1)*nx-i-2 -1,1)
    mp1=np.delete(mp1,3*(ny-2)*nx+1-i+2 -1,0)
    mp1=np.delete(mp1,3*(ny-2)*nx+1-i+2 -1,1)
    mp1=np.delete(mp1,3*(ny-2)*nx+1-i+1 -1,0)
    mp1=np.delete(mp1,3*(ny-2)*nx+1-i+1 -1,1)
    mp1=np.delete(mp1,3*(ny-2)*nx+1-i -1,0)
    mp1=np.delete(mp1,3*(ny-2)*nx+1-i -1,1)
for i in range(0,3*nx-1 +1):
    mp1=np.delete(mp1,3*nx-i -1,0)
    mp1=np.delete(mp1,3*nx-i -1,1)

# compute eigenfrequencies
eival=np.linalg.eigvals(np.dot(np.linalg.inv(mp1),kp1))
eigf=(np.sqrt(Di/(rost*h))/(2*np.pi))*np.sqrt(eival)
print(np.sort(eigf)[0:12])
