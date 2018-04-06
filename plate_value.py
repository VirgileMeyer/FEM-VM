import numpy as np
# Medium values
c=343.                  # Sound celerity (in m/s)
ro=1.2                  # Air density (in kg)

# Geometrical values
lx=0.3                 # Plate dimension along X (in m)
ly=0.6                  # Plate dimension along Y (in m)

nx=int(35)               # Number of nodes along X
ny=int(65)              # Number of nodes along Y

dx=lx/(nx-1)
dy=ly/(ny-1)



freq=100.                # Frequency (in Hz)


# Material values

# Aluminium
Ei=70.*10**9             # Young modulud (in Pa)
nu=0.3                   # Poisson coefficient
rost=2700.               # Surface density (in kg/m**2)
h=0.0001                 # Plate thickness (in m)

# Steel
##rost=7850
##Ei=210.*10**9
##nu=0.3
##h=0.0001

Di=Ei*h**3/(12*(1-nu**2))


#elementary matrices
aa=dx/2
bb=dy/2

#kep=np.zeros(12)
##kep=np.full((12,12),1)

kep=np.array([[(8*aa*bb*(10*dx**4 + 10*dy**4 + dx**2*dy**2*(7 - 2*nu)))/(5.*dx**4*dy**4), \
(4*aa*bb*(10*dx**2 + dy**2*(1 + 4*nu)))/(5.*dx**2*dy**3),(-4*aa*bb*(10*dy**2 + dx**2*(1 + 4*nu)))/(5.*dx**3*dy**2), \
(8*aa*bb*(5*dx**4 - 10*dy**4 + dx**2*dy**2*(-7 + 2*nu)))/(5.*dx**4*dy**4),(-4*aa*bb*(-5*dx**2 + dy**2*(1 + 4*nu)))/(5.*dx**2*dy**3), \
(aa*bb*(-40*dy**2 + 4*dx**2*(-1 + nu)))/(5.*dx**3*dy**2),(-8*aa*bb*(5*dx**4 + 5*dy**4 + dx**2*dy**2*(-7 + 2*nu)))/(5.*dx**4*dy**4), \
(4*aa*bb*(5*dx**2 + dy**2*(-1 + nu)))/(5.*dx**2*dy**3),(-4*aa*bb*(5*dy**2 + dx**2*(-1 + nu)))/(5.*dx**3*dy**2), \
(-8*aa*bb*(10*dx**4 - 5*dy**4 + dx**2*dy**2*(7 - 2*nu)))/(5.*dx**4*dy**4),(4*aa*bb*(10*dx**2 - dy**2*(-1 + nu)))/(5.*dx**2*dy**3), \
(4*aa*bb*(-5*dy**2 + dx**2*(1 + 4*nu)))/(5.*dx**3*dy**2)], \
[(4*aa*bb*(10*dx**2 + dy**2*(1 + 4*nu)))/(5.*dx**2*dy**3),(16*aa*bb*(5*dx**2 - dy**2*(-1 + nu)))/(15.*dx**2*dy**2), \
(-4*aa*bb*nu)/(dx*dy),(-4*aa*bb*(-5*dx**2 + dy**2*(1 + 4*nu)))/(5.*dx**2*dy**3), \
(8*aa*bb*(5*dx**2 + 2*dy**2*(-1 + nu)))/(15.*dx**2*dy**2),0,(-4*aa*bb*(5*dx**2 + dy**2*(-1 + nu)))/(5.*dx**2*dy**3), \
(4*aa*bb*(5/dy**2 + (1 - nu)/dx**2))/15.,0,(-4*aa*bb*(10*dx**2 - dy**2*(-1 + nu)))/(5.*dx**2*dy**3), \
(4*aa*bb*(10/dy**2 + (-1 + nu)/dx**2))/15.,0],[(-4*aa*bb*(10*dy**2 + dx**2*(1 + 4*nu)))/(5.*dx**3*dy**2),(-4*aa*bb*nu)/(dx*dy), \
(16*aa*bb*(5*dy**2 - dx**2*(-1 + nu)))/(15.*dx**2*dy**2),(aa*bb*(40*dy**2 - 4*dx**2*(-1 + nu)))/(5.*dx**3*dy**2),0, \
(4*aa*bb*(10/dx**2 + (-1 + nu)/dy**2))/15.,(4*aa*bb*(5*dy**2 + dx**2*(-1 + nu)))/(5.*dx**3*dy**2),0, \
(4*aa*bb*(5/dx**2 + (1 - nu)/dy**2))/15.,(4*aa*bb*(-5*dy**2 + dx**2*(1 + 4*nu)))/(5.*dx**3*dy**2),0, \
(8*aa*bb*(5*dy**2 + 2*dx**2*(-1 + nu)))/(15.*dx**2*dy**2)], \
[(8*aa*bb*(5*dx**4 - 10*dy**4 + dx**2*dy**2*(-7 + 2*nu)))/(5.*dx**4*dy**4),(-4*aa*bb*(-5*dx**2 + dy**2*(1 + 4*nu)))/(5.*dx**2*dy**3), \
(aa*bb*(40*dy**2 - 4*dx**2*(-1 + nu)))/(5.*dx**3*dy**2),(8*aa*bb*(10*dx**4 + 10*dy**4 + dx**2*dy**2*(7 - 2*nu)))/(5.*dx**4*dy**4), \
(4*aa*bb*(10*dx**2 + dy**2*(1 + 4*nu)))/(5.*dx**2*dy**3),(4*aa*bb*(10*dy**2 + dx**2*(1 + 4*nu)))/(5.*dx**3*dy**2), \
(-8*aa*bb*(10*dx**4 - 5*dy**4 + dx**2*dy**2*(7 - 2*nu)))/(5.*dx**4*dy**4),(4*aa*bb*(10*dx**2 - dy**2*(-1 + nu)))/(5.*dx**2*dy**3), \
(-4*aa*bb*(-5*dy**2 + dx**2*(1 + 4*nu)))/(5.*dx**3*dy**2),(-8*aa*bb*(5*dx**4 + 5*dy**4 + dx**2*dy**2*(-7 + 2*nu)))/(5.*dx**4*dy**4), \
(4*aa*bb*(5*dx**2 + dy**2*(-1 + nu)))/(5.*dx**2*dy**3),(4*aa*bb*(5*dy**2 + dx**2*(-1 + nu)))/(5.*dx**3*dy**2)], \
[(-4*aa*bb*(-5*dx**2 + dy**2*(1 + 4*nu)))/(5.*dx**2*dy**3),(8*aa*bb*(5*dx**2 + 2*dy**2*(-1 + nu)))/(15.*dx**2*dy**2),0, \
(4*aa*bb*(10*dx**2 + dy**2*(1 + 4*nu)))/(5.*dx**2*dy**3),(16*aa*bb*(5*dx**2 - dy**2*(-1 + nu)))/(15.*dx**2*dy**2),(4*aa*bb*nu)/(dx*dy), \
(-4*aa*bb*(10*dx**2 - dy**2*(-1 + nu)))/(5.*dx**2*dy**3),(4*aa*bb*(10/dy**2 + (-1 + nu)/dx**2))/15.,0, \
(-4*aa*bb*(5*dx**2 + dy**2*(-1 + nu)))/(5.*dx**2*dy**3),(4*aa*bb*(5/dy**2 + (1 - nu)/dx**2))/15.,0], \
[(aa*bb*(-40*dy**2 + 4*dx**2*(-1 + nu)))/(5.*dx**3*dy**2),0,(4*aa*bb*(10/dx**2 + (-1 + nu)/dy**2))/15., \
(4*aa*bb*(10*dy**2 + dx**2*(1 + 4*nu)))/(5.*dx**3*dy**2),(4*aa*bb*nu)/(dx*dy),(16*aa*bb*(5*dy**2 - dx**2*(-1 + nu)))/(15.*dx**2*dy**2), \
(-4*aa*bb*(-5*dy**2 + dx**2*(1 + 4*nu)))/(5.*dx**3*dy**2),0,(8*aa*bb*(5*dy**2 + 2*dx**2*(-1 + nu)))/(15.*dx**2*dy**2), \
(-4*aa*bb*(5*dy**2 + dx**2*(-1 + nu)))/(5.*dx**3*dy**2),0,(4*aa*bb*(5/dx**2 + (1 - nu)/dy**2))/15.], \
[(-8*aa*bb*(5*dx**4 + 5*dy**4 + dx**2*dy**2*(-7 + 2*nu)))/(5.*dx**4*dy**4),(-4*aa*bb*(5*dx**2 + dy**2*(-1 + nu)))/(5.*dx**2*dy**3), \
(4*aa*bb*(5*dy**2 + dx**2*(-1 + nu)))/(5.*dx**3*dy**2),(-8*aa*bb*(10*dx**4 - 5*dy**4 + dx**2*dy**2*(7 - 2*nu)))/(5.*dx**4*dy**4), \
(-4*aa*bb*(10*dx**2 - dy**2*(-1 + nu)))/(5.*dx**2*dy**3),(-4*aa*bb*(-5*dy**2 + dx**2*(1 + 4*nu)))/(5.*dx**3*dy**2), \
(8*aa*bb*(10*dx**4 + 10*dy**4 + dx**2*dy**2*(7 - 2*nu)))/(5.*dx**4*dy**4),(-4*aa*bb*(10*dx**2 + dy**2*(1 + 4*nu)))/(5.*dx**2*dy**3), \
(4*aa*bb*(10*dy**2 + dx**2*(1 + 4*nu)))/(5.*dx**3*dy**2),(8*aa*bb*(5*dx**4 - 10*dy**4 + dx**2*dy**2*(-7 + 2*nu)))/(5.*dx**4*dy**4), \
(4*aa*bb*(-5*dx**2 + dy**2*(1 + 4*nu)))/(5.*dx**2*dy**3),(aa*bb*(40*dy**2 - 4*dx**2*(-1 + nu)))/(5.*dx**3*dy**2)], \
[(4*aa*bb*(5*dx**2 + dy**2*(-1 + nu)))/(5.*dx**2*dy**3),(4*aa*bb*(5/dy**2 + (1 - nu)/dx**2))/15.,0, \
(4*aa*bb*(10*dx**2 - dy**2*(-1 + nu)))/(5.*dx**2*dy**3),(4*aa*bb*(10/dy**2 + (-1 + nu)/dx**2))/15.,0, \
(-4*aa*bb*(10*dx**2 + dy**2*(1 + 4*nu)))/(5.*dx**2*dy**3),(16*aa*bb*(5*dx**2 - dy**2*(-1 + nu)))/(15.*dx**2*dy**2), \
(-4*aa*bb*nu)/(dx*dy),(4*aa*bb*(-5*dx**2 + dy**2*(1 + 4*nu)))/(5.*dx**2*dy**3), \
(8*aa*bb*(5*dx**2 + 2*dy**2*(-1 + nu)))/(15.*dx**2*dy**2),0], \
[(-4*aa*bb*(5*dy**2 + dx**2*(-1 + nu)))/(5.*dx**3*dy**2),0,(4*aa*bb*(5/dx**2 + (1 - nu)/dy**2))/15., \
(-4*aa*bb*(-5*dy**2 + dx**2*(1 + 4*nu)))/(5.*dx**3*dy**2),0,(8*aa*bb*(5*dy**2 + 2*dx**2*(-1 + nu)))/(15.*dx**2*dy**2), \
(4*aa*bb*(10*dy**2 + dx**2*(1 + 4*nu)))/(5.*dx**3*dy**2),(-4*aa*bb*nu)/(dx*dy),(16*aa*bb*(5*dy**2 - dx**2*(-1 + nu)))/(15.*dx**2*dy**2), \
(aa*bb*(-40*dy**2 + 4*dx**2*(-1 + nu)))/(5.*dx**3*dy**2),0,(4*aa*bb*(10/dx**2 + (-1 + nu)/dy**2))/15.], \
[(-8*aa*bb*(10*dx**4 - 5*dy**4 + dx**2*dy**2*(7 - 2*nu)))/(5.*dx**4*dy**4),(-4*aa*bb*(10*dx**2 - dy**2*(-1 + nu)))/(5.*dx**2*dy**3), \
(4*aa*bb*(-5*dy**2 + dx**2*(1 + 4*nu)))/(5.*dx**3*dy**2),(-8*aa*bb*(5*dx**4 + 5*dy**4 + dx**2*dy**2*(-7 + 2*nu)))/(5.*dx**4*dy**4), \
(-4*aa*bb*(5*dx**2 + dy**2*(-1 + nu)))/(5.*dx**2*dy**3),(-4*aa*bb*(5*dy**2 + dx**2*(-1 + nu)))/(5.*dx**3*dy**2), \
(8*aa*bb*(5*dx**4 - 10*dy**4 + dx**2*dy**2*(-7 + 2*nu)))/(5.*dx**4*dy**4),(4*aa*bb*(-5*dx**2 + dy**2*(1 + 4*nu)))/(5.*dx**2*dy**3), \
(aa*bb*(-40*dy**2 + 4*dx**2*(-1 + nu)))/(5.*dx**3*dy**2),(8*aa*bb*(10*dx**4 + 10*dy**4 + dx**2*dy**2*(7 - 2*nu)))/(5.*dx**4*dy**4), \
(-4*aa*bb*(10*dx**2 + dy**2*(1 + 4*nu)))/(5.*dx**2*dy**3),(-4*aa*bb*(10*dy**2 + dx**2*(1 + 4*nu)))/(5.*dx**3*dy**2)], \
[(4*aa*bb*(10*dx**2 - dy**2*(-1 + nu)))/(5.*dx**2*dy**3),(4*aa*bb*(10/dy**2 + (-1 + nu)/dx**2))/15.,0, \
(4*aa*bb*(5*dx**2 + dy**2*(-1 + nu)))/(5.*dx**2*dy**3),(4*aa*bb*(5/dy**2 + (1 - nu)/dx**2))/15.,0, \
(4*aa*bb*(-5*dx**2 + dy**2*(1 + 4*nu)))/(5.*dx**2*dy**3),(8*aa*bb*(5*dx**2 + 2*dy**2*(-1 + nu)))/(15.*dx**2*dy**2),0, \
(-4*aa*bb*(10*dx**2 + dy**2*(1 + 4*nu)))/(5.*dx**2*dy**3),(16*aa*bb*(5*dx**2 - dy**2*(-1 + nu)))/(15.*dx**2*dy**2),(4*aa*bb*nu)/(dx*dy)] \
,[(4*aa*bb*(-5*dy**2 + dx**2*(1 + 4*nu)))/(5.*dx**3*dy**2),0,(8*aa*bb*(5*dy**2 + 2*dx**2*(-1 + nu)))/(15.*dx**2*dy**2), \
(4*aa*bb*(5*dy**2 + dx**2*(-1 + nu)))/(5.*dx**3*dy**2),0,(4*aa*bb*(5/dx**2 + (1 - nu)/dy**2))/15., \
(aa*bb*(40*dy**2 - 4*dx**2*(-1 + nu)))/(5.*dx**3*dy**2),0,(4*aa*bb*(10/dx**2 + (-1 + nu)/dy**2))/15., \
(-4*aa*bb*(10*dy**2 + dx**2*(1 + 4*nu)))/(5.*dx**3*dy**2),(4*aa*bb*nu)/(dx*dy),(16*aa*bb*(5*dy**2 - dx**2*(-1 + nu)))/(15.*dx**2*dy**2)]])

mep=np.array([[(1727*aa*bb)/3150.,(461*aa*bb**2)/3150., \
(-461*aa**2*bb)/3150.,(613*aa*bb)/3150.,(199*aa*bb**2)/3150., \
(137*aa**2*bb)/1575.,(197*aa*bb)/3150.,(-58*aa*bb**2)/1575., \
(58*aa**2*bb)/1575.,(613*aa*bb)/3150.,(-137*aa*bb**2)/1575., \
(-199*aa**2*bb)/3150.], \
[(461*aa*bb**2)/3150.,(16*aa*bb**3)/315.,-(aa**2*bb**2)/25., \
(199*aa*bb**2)/3150.,(8*aa*bb**3)/315.,(2*aa**2*bb**2)/75., \
(58*aa*bb**2)/1575.,(-2*aa*bb**3)/105.,(4*aa**2*bb**2)/225., \
(137*aa*bb**2)/1575.,(-4*aa*bb**3)/105.,(-2*aa**2*bb**2)/75.], \
[(-461*aa**2*bb)/3150.,-(aa**2*bb**2)/25., \
(16*aa**3*bb)/315.,(-137*aa**2*bb)/1575.,(-2*aa**2*bb**2)/75., \
(-4*aa**3*bb)/105.,(-58*aa**2*bb)/1575.,(4*aa**2*bb**2)/225., \
(-2*aa**3*bb)/105.,(-199*aa**2*bb)/3150.,(2*aa**2*bb**2)/75., \
(8*aa**3*bb)/315.],[(613*aa*bb)/3150.,(199*aa*bb**2)/3150., \
(-137*aa**2*bb)/1575.,(1727*aa*bb)/3150.,(461*aa*bb**2)/3150., \
(461*aa**2*bb)/3150.,(613*aa*bb)/3150.,(-137*aa*bb**2)/1575., \
(199*aa**2*bb)/3150.,(197*aa*bb)/3150.,(-58*aa*bb**2)/1575., \
(-58*aa**2*bb)/1575.], \
[(199*aa*bb**2)/3150.,(8*aa*bb**3)/315., \
(-2*aa**2*bb**2)/75.,(461*aa*bb**2)/3150.,(16*aa*bb**3)/315., \
(aa**2*bb**2)/25.,(137*aa*bb**2)/1575.,(-4*aa*bb**3)/105., \
(2*aa**2*bb**2)/75.,(58*aa*bb**2)/1575.,(-2*aa*bb**3)/105., \
(-4*aa**2*bb**2)/225.], \
[(137*aa**2*bb)/1575.,(2*aa**2*bb**2)/75., \
(-4*aa**3*bb)/105.,(461*aa**2*bb)/3150.,(aa**2*bb**2)/25., \
(16*aa**3*bb)/315.,(199*aa**2*bb)/3150.,(-2*aa**2*bb**2)/75., \
(8*aa**3*bb)/315.,(58*aa**2*bb)/1575.,(-4*aa**2*bb**2)/225., \
(-2*aa**3*bb)/105.], \
[(197*aa*bb)/3150.,(58*aa*bb**2)/1575.,(-58*aa**2*bb)/1575., \
(613*aa*bb)/3150.,(137*aa*bb**2)/1575.,(199*aa**2*bb)/3150., \
(1727*aa*bb)/3150.,(-461*aa*bb**2)/3150.,(461*aa**2*bb)/3150., \
(613*aa*bb)/3150.,(-199*aa*bb**2)/3150.,(-137*aa**2*bb)/1575.], \
[(-58*aa*bb**2)/1575.,(-2*aa*bb**3)/105., \
(4*aa**2*bb**2)/225.,(-137*aa*bb**2)/1575.,(-4*aa*bb**3)/105., \
(-2*aa**2*bb**2)/75.,(-461*aa*bb**2)/3150.,(16*aa*bb**3)/315., \
-(aa**2*bb**2)/25.,(-199*aa*bb**2)/3150.,(8*aa*bb**3)/315., \
(2*aa**2*bb**2)/75.], \
[(58*aa**2*bb)/1575.,(4*aa**2*bb**2)/225., \
(-2*aa**3*bb)/105.,(199*aa**2*bb)/3150.,(2*aa**2*bb**2)/75., \
(8*aa**3*bb)/315.,(461*aa**2*bb)/3150.,-(aa**2*bb**2)/25., \
(16*aa**3*bb)/315.,(137*aa**2*bb)/1575.,(-2*aa**2*bb**2)/75., \
(-4*aa**3*bb)/105.], \
[(613*aa*bb)/3150.,(137*aa*bb**2)/1575., \
(-199*aa**2*bb)/3150.,(197*aa*bb)/3150.,(58*aa*bb**2)/1575., \
(58*aa**2*bb)/1575.,(613*aa*bb)/3150.,(-199*aa*bb**2)/3150., \
(137*aa**2*bb)/1575.,(1727*aa*bb)/3150.,(-461*aa*bb**2)/3150., \
(-461*aa**2*bb)/3150.], \
[(-137*aa*bb**2)/1575.,(-4*aa*bb**3)/105., \
(2*aa**2*bb**2)/75.,(-58*aa*bb**2)/1575.,(-2*aa*bb**3)/105., \
(-4*aa**2*bb**2)/225.,(-199*aa*bb**2)/3150.,(8*aa*bb**3)/315., \
(-2*aa**2*bb**2)/75.,(-461*aa*bb**2)/3150.,(16*aa*bb**3)/315., \
(aa**2*bb**2)/25.],[(-199*aa**2*bb)/3150., \
(-2*aa**2*bb**2)/75.,(8*aa**3*bb)/315.,(-58*aa**2*bb)/1575., \
(-4*aa**2*bb**2)/225.,(-2*aa**3*bb)/105.,(-137*aa**2*bb)/1575., \
(2*aa**2*bb**2)/75.,(-4*aa**3*bb)/105.,(-461*aa**2*bb)/3150., \
(aa**2*bb**2)/25.,(16*aa**3*bb)/315.]])
