“ !!!! work in progress !!!! ”
“ I am currently working on the validation using a PDE solver “
“ The eigenvalues seem to begin good agreement with semi-empirical method found in ESDU 75030 ”




This Code written in Python will compute the displacement from 
the vibration of a rectangular Kirchhoff–Love plate using the finite element method

To use it, it requires the libraries numpy and matplotlib

Boundary Condition : 
 - ‘Clamped *4’

The elements used are of the “plate type” : 4 nodes with 3 degrees of freedom at each node
(1 for displacement, 2 for in-plane rotation)


Input values (e.g.. Young’Modulus, geometry, number of nodes…) are written in ‘plate_value.py’
(beginning of the file)

Run the file ‘Plate_EF.py’ to plot the displacement field

Run the file ‘eigenfreq.py’ to compute the eigenfrequencies (Natural frequencies)