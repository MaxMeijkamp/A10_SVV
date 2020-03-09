### Importing required packages
import numpy as np
import Energy
import Stiffness
import Stress
import math as m

######################## Part I - parameters as in assignment #######################################
aircraft = "CRJ700" # Write either A320, F100, CRJ700 or Do228 (bear in mind capitals); this is used for aerodynamic loading
Ca = 0.484  # m
la = 1.691  # m
x1 = 0.149  # m
x2 = 0.554  # m
x3 = 1.541  # m
xa = 0.272  # m
ha = 0.173  # m
tsk = 0.0011  # m
tsp = 0.0025  # m
tst = 0.0012  # m
hst = 0.014   # m
wst = 0.018   # m
nst = 13  # -
d1 = 0.00681  # m
d3 = 0.02030  # m
theta = m.radians(26)  # rad
P = 37.9*1000  # N

######################## Part II - bending stiffness calculations #######################################
### Create the cross-section object
"""" Note that only the cross-sectional geometry is put through. Furthermore, you cannot disable this line, as the 
cross-section object is used in subsequent calculations"""
crosssection = Stiffness.Crosssection(nst,Ca,ha,tsk,tsp,tst,hst,wst)  # Define the cross-section

### Primary functions
""" If you desire, you may disable this line, and manually overwrite the values listed between lines 45-50"""
crosssection.compute_bending_properties()   # Run the calculations

### Auxiliary functions
""" A plot of the cross-section, to inspect that the stringers have been placed correctly, and 
that the position of the centroid makes sense. """
#crosssection.plot_crosssection()     # Plot the cross-section; blue cross is the centroid, red crosses are stringers

### Access to important results
"""" If you desire, you can manually overwrite these values. """
stcoord = crosssection.stcoord            # array containing stringer coordinates
_ = crosssection.totarea            # total cross-section area
yc = crosssection.yc                 # y-coordinate of the centroid
zc = crosssection.zc                 # z-coordinate of the centroid
Iyy = crosssection.Iyy                # moment of inertia about y-axis
Izz = crosssection.Izz                # moment of inertia about z-axis
#print(zc)

######################## Part III - Torsional stiffness calculations #######################################
### Primary functions
""" If you desire, you may disable this line, and manually overwrite the values listed between lines 60-62"""
crosssection.compute_shearcenter()   # Run the calculations
crosssection.compute_torsionalstiffness()   # Run the calculations

h = crosssection.ha / 2.
A1 = m.pi * h ** 2 / 2.
A2 = (crosssection.Ca - h) * h

A = np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])
b = np.array([0., 0., 0.])

### First row
A[0, 0] = 2. * A1
A[0, 1] = 2. * A2
b[0] = 1

### Second row
A[1, 0] = (h * m.pi / crosssection.tsk + 2 * h / crosssection.tsp) / (2 * A1)
A[1, 1] = (-2 * h / crosssection.tsp) / (2 * A1)
A[1, 2] = -1.
b[1] = 0.

### Third row
A[2, 0] = (-2 * h / crosssection.tsp) / (2 * A2)
A[2, 1] = (2 * crosssection.lsk / crosssection.tsk + 2 * h / crosssection.tsp) / (2 * A2)
A[2, 2] = -1
b[2] = 0.

solution = np.linalg.solve(A, b)
crosssection.J = 1. / solution[-1]

### Access to important results
"""" If you desire, you can manually overwrite these values. """
ysc = crosssection.ysc                 # y-coordinate of the centroid
zsc = crosssection.zsc                 # z-coordinate of the centroid
#print("ysc, zsc = ",ysc, zsc+h)

_ = crosssection.J                   # torsional constant

######################## Part IV - Deflection calculations #######################################
### Definition of additional parameters
N = 20     # Number of basis functions to use in Rayleigh-Ritz method (total number of coefficients is 3*N)
E = 73.1*10**9       # E-modulus (Pa)
G = 28*10**9       # G-modulus (Pa)

### Create the aileron object
""" Merges the cross-sectional properties with the spanwise properties (length and material properties)"""
aileron = Energy.Beam(la,crosssection,N,E,G)

"""Define your boundary conditions; see manual for explanations. The shown boundary conditions are the boundary
conditions for the aileron as described in the assignment."""
aileron.addbcss(x1,0.,-ha/2.,-theta,d1)
aileron.addbcss(x1,0.,-ha/2.,m.pi/2-theta,0)
aileron.addbcss(x2,0.,-ha/2.,0,0)
aileron.addbcss(x2,0.,-ha/2.,m.pi/2,0.)
aileron.addbcss(x3,0.,-ha/2.,-theta,d3)
aileron.addbcss(x3,0.,-ha/2.,m.pi/2-theta,0)
aileron.addbcss(x2-xa/2.,ha/2.,0,m.pi/2.-theta,0)

""""Define your applied loading; see manual for explanations."""
aileron.addfpl(x2+xa/2.,ha/2.,0,m.pi/2.-theta,-P)

### Primary functions
""" The following line computes the deflections. If you do not want to include the aerodynamic loading, simply write
aileron.compute_deflections(). Note that the aerodynamic loading significantly slows down the program.
If you do want to include the aerodynamic loading, let the variable aircraft (see part I) equal "A320", "F100", "CRJ700",
"Do228", whichever you want to include. 
Note that the name should be spelled exactly as listed above. Note that if the aircraft you write is inconsistent with the
geometry you define at the beginning of this file, the program will not return an error, but will simply produce bogus
results."""
aileron.compute_deflections() ### Switch aerodynamic loading to the aircraft that is being considered

### Auxiliary functions
"""" A number of auxiliary functions and results are given to you. """

## Simplistic plotting procedures for a first check
#aileron.plotv()             # Plot the deflections in y-direction, its derivative, the bending moment about the z-axis, and the shear force in y.
#aileron.plotw()             # Plot the deflections in z-direction, its derivative, the bending moment about the y-axis, and the shear force in z.
#aileron.plotphi()           # Plot the twist distribution, the torque and the distributed torque.

## For custom post-processing of the solution
x = np.linspace(0,la,num = 100)  # Subsequent functions accept numpy-arrays
# Compute the deflections
vdef, wdef, phidef = aileron.eval(x)       # Compute the three deflections
_, _, _ = aileron.fdeval(x)     # Compute their their first order derivative
_, _, _ = aileron.sdeval(x)     # Compute their their second order derivative
_, _, _ = aileron.tdeval(x)     # Compute their their third order derivative
#print(x, vdef, wdef, phidef)

# Compute the loading
_ = aileron.Sy(x)               # Compute the shear force in y
_ = aileron.Sz(x)               # Compute the shear force in z
_ = aileron.My(x)               # Compute the moment around the y-axis
_ = aileron.Mz(x)               # Compute the moment around the z-axis
_ = aileron.T(x)                # Compute the torque
_ = aileron.tau(x)              # Compute the distributed torque

## Value of the total potential energy
_ = aileron.cPI()               # Compute the total potential energy of the beam for the computed solution.
i,k = 0, 1.01                   # Parameters for next line
_ = aileron.cPI_vary_coef(i,k)  # Multiply the i-th coefficient of hat(alpha) by a factor k and compute the corresponding TPE.


### Access to important results
_ = aileron.Na              # Number of coefficients used to approximate v(x)
_ = aileron.Nb              # Number of coefficients used to approximate w(x)
_ = aileron.Nc              # Number of coefficients used to approximate phi(x)
_ = aileron.nbc             # Total number of boundary conditions
_ = aileron.nbcv            # Number of coefficients used for boundary conditions for v; equal to N - Na
_ = aileron.nbcw            # Number of coefficients used for boundary conditions for w; equal to N - Nb
_ = aileron.nbct            # Number of coefficients used for boundary conditions for phi; equal to N - Nc

_ = aileron.Ha              # H_a matrix
_ = aileron.Hb              # H_b matrix
_ = aileron.Hc              # H_c matrix

_ = aileron.Ua              # Upsilon_a matrix
_ = aileron.Ub              # Upsilon_b matrix
_ = aileron.Uc              # Upsilon_c matrix

_ = aileron.K1              # K_{1,a}/K_{1,b} matrix
_ = aileron.C1              # K_{1,c} matrix
_ = aileron.K2a             # K_{2,a} matrix
_ = aileron.K2b             # K_{2,b} matrix
_ = aileron.C2              # K_{2,c} vector
_ = aileron.F               # F vector

_ = aileron.LHS             # Left-hand-side matrix
_ = aileron.RHS             # Right-hand-side vector

_ = aileron.sol.coef        # Resulting coefficients, collected in bar(alpha) (and thus include both the 'boundary' and 'external' coefficients.

######################## Part V - Stress calculations #######################################
### Create the stress state object, which will contain information about the stresses of the aileron
Stressobject = Stress.Stressstate(crosssection)

### Define the forces and moments for which you want to know the stress distributions
xmax = []
ymax = []
zmax = []
qmax = []
qmin = []
sssmax = []
sssmin = []
vmmax = []
vmmin = []
verstress = np.array([])

for x in range(0,5):
    verstresstemp = []
    xtemp = x/1000
    Sy = aileron.Sy(xtemp)
    Sz = aileron.Sz(xtemp)
    My = aileron.My(xtemp)
    Mz = aileron.Mz(xtemp)
    T = aileron.T(xtemp)

    ### Primary functions
    """"The following line should never be disabled, as its results are used in the auxiliary functions"""
    Stressobject.compute_unitstressdistributions()

    h = Stressobject.ha / 2.
    A1 = m.pi * h ** 2 / 2.
    A2 = (Stressobject.Ca - h) * h

    A = np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])
    b = np.array([0., 0., 0.])

    ### First row
    A[0, 0] = 2. * A1
    A[0, 1] = 2. * A2
    b[0] = -1

    ### Second row
    A[1, 0] = (h * m.pi / Stressobject.tsk + 2 * h / Stressobject.tsp) / (2 * A1)
    A[1, 1] = (-2 * h / Stressobject.tsp) / (2 * A1)
    A[1, 2] = -1.
    b[1] = 0.

    ### Third row
    A[2, 0] = (-2 * h / Stressobject.tsp) / (2 * A2)
    A[2, 1] = (2 * Stressobject.lsk / Stressobject.tsk + 2 * h / Stressobject.tsp) / (2 * A2)
    A[2, 2] = -1
    b[2] = 0.

    solution = np.linalg.solve(A, b)
    Stressobject.Tq1f = solution[0]
    Stressobject.Tq2f = -solution[0] + solution[1]
    Stressobject.Tq3f = solution[1]
    Stressobject.Tq4f = solution[1]
    Stressobject.Tq5f = -solution[0] + solution[1]
    Stressobject.Tq6f = solution[0]

    ### Auxiliary functions
    Stressobject.compute_stressdistributions(Sy,Sz,My,Mz,T)

    ### Some plotting functions
    #Stressobject.plot_shearflowdistributions()
    #Stressobject.plot_directstressdistributions()
    #Stressobject.plot_vonmisesstressdistributions()

    ### Access to important results
    theta = np.linspace(0,m.pi/2,num = 100)
    q1 = Stressobject.q1f(theta)             # Compute the shear flow distribution in region 1
    sss1 = Stressobject.sigma1f(theta)         # Compute the direct stress distribution in region 1
    vm1 = Stressobject.vm1(theta)             # Compute the Von Mises stress distribution in region 1
    z1, y1 = Stressobject.coord1(theta)       # Compute the z,y-coordinates for region 1
    verstresstemp = np.array([xtemp,y1,z1,q1,sss1,vm1])
    verstress.append(verstresstemp)

    y = np.linspace(0,ha/2.,num = 100)
    q2 = Stressobject.q2f(y)             # Compute the shear flow distribution in region 3
    sss2 = Stressobject.sigma2f(y)         # Compute the direct stress distribution in region 3
    vm2 = Stressobject.vm2(y)             # Compute the Von Mises stress distribution in region 3
    z2, y2 = Stressobject.coord2(y)       # Compute the z,y-coordinates for region 3
    verstresstemp = [xtemp,y2,z2,q2,sss2,vm2]
    verstress.append(verstresstemp)

    s = np.linspace(0,m.sqrt((Ca-ha/2.)**2+(ha/2.)**2),num = 100)
    q3 = Stressobject.q3f(s)             # Compute the shear flow distribution in region 4
    sss3 = Stressobject.sigma3f(s)         # Compute the direct stress distribution in region 4
    vm3 = Stressobject.vm3(s)             # Compute the Von Mises stress distribution in region 4
    z3, y3 = Stressobject.coord3(s)       # Compute the z,y-coordinates for region 4
    verstresstemp = [xtemp,y3,z3,q3,sss3,vm3]
    verstress.append(verstresstemp)

    s = np.linspace(0,m.sqrt((Ca-ha/2.)**2+(ha/2.)**2),num = 100)
    q4 = Stressobject.q4f(s)             # Compute the shear flow distribution in region 4
    sss4 = Stressobject.sigma4f(s)         # Compute the direct stress distribution in region 4
    vm4 = Stressobject.vm4(s)             # Compute the Von Mises stress distribution in region 4
    z4, y4 = Stressobject.coord4(s)       # Compute the z,y-coordinates for region 4
    verstresstemp = [xtemp,y4,z4,q4,sss4,vm4]
    verstress.append(verstresstemp)

    y = np.linspace(0,ha/2.,num = 100)
    q5 = Stressobject.q5f(y)             # Compute the shear flow distribution in region 5
    sss5 = Stressobject.sigma5f(y)         # Compute the direct stress distribution in region 5
    vm5 = Stressobject.vm5(y)             # Compute the Von Mises stress distribution in region 5
    z5, y5 = Stressobject.coord5(y)       # Compute the z,y-coordinates for region 5
    verstresstemp = [xtemp,y5,z5,q5,sss5,vm5]
    verstress.append(verstresstemp)

    theta = np.linspace(-m.pi/2,0,num = 100)
    q6 = Stressobject.q6f(theta)             # Compute the shear flow distribution in region 6
    sss6 = Stressobject.sigma6f(theta)         # Compute the direct stress distribution in region 6
    vm6 = Stressobject.vm6(theta)             # Compute the Von Mises stress distribution in region 6
    z6, y6 = Stressobject.coord6(theta)       # Compute the z,y-coordinates for region 6
    verstresstemp = [xtemp,y6,z6,q6,sss6,vm6]
    verstress.append(verstresstemp)

    #print(xtemp, max(vm1), min(vm1))
    #print(xtemp, max(vm2), min(vm2))
    #print(xtemp, max(vm3), min(vm3))
    #print(xtemp, max(vm4), min(vm4))
    #print(xtemp, max(vm5), min(vm5))
    #print(xtemp, max(vm6), min(vm6))
    qmaxtemp = max(max(q1), max(q2), max(q3), max(q4), max(q5), max(q6))
    qmintemp = min(min(q1), min(q2), min(q3), min(q4), min(q5), min(q6))
    sssmaxtemp = max(max(sss1), max(sss2), max(sss3), max(sss4), max(sss5), max(sss6))
    sssmintemp = min(min(sss1), min(sss2), min(sss3), min(sss4), min(sss5), min(sss6))
    vmmaxtemp = max(max(vm1), max(vm2), max(vm3), max(vm4), max(vm5), max(vm6))
    vmmintemp = min(min(vm1), min(vm2), min(vm3), min(vm4), min(vm5), min(vm6))

    xmax.append(xtemp)
    qmax.append(qmaxtemp)
    qmin.append(qmintemp)
    sssmax.append(sssmaxtemp)
    sssmin.append(sssmintemp)

    vmmax.append(vmmaxtemp)
    vmmin.append(vmmintemp)

i = qmax.index(max(qmax))
j = qmin.index(min(qmin))
k = sssmax.index(max(sssmax))
l = sssmin.index(min(sssmin))
r = vmmax.index(max(vmmax))
n = vmmin.index(min(vmmin))

print(xmax[j],qmin[j],xmax[i],qmax[i])
print(xmax[l],sssmin[l],xmax[k],sssmax[k])
print(xmax[n],vmmin[n],xmax[r],vmmax[r])

xstress = [xmax[i],xmax[k],xmax[l],xmax[n]]
print(xstress)
for p in range(4):
    print("wollah", p)
    x = xstress[i]
    Sy = aileron.Sy(xtemp)
    Sz = aileron.Sz(xtemp)
    My = aileron.My(xtemp)
    Mz = aileron.Mz(xtemp)
    T = aileron.T(xtemp)
    Stressobject.compute_unitstressdistributions()
    h = Stressobject.ha / 2.
    A1 = m.pi * h ** 2 / 2.
    A2 = (Stressobject.Ca - h) * h

    A = np.array([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])
    b = np.array([0., 0., 0.])

    ### First row
    A[0, 0] = 2. * A1
    A[0, 1] = 2. * A2
    b[0] = -1

    ### Second row
    A[1, 0] = (h * m.pi / Stressobject.tsk + 2 * h / Stressobject.tsp) / (2 * A1)
    A[1, 1] = (-2 * h / Stressobject.tsp) / (2 * A1)
    A[1, 2] = -1.
    b[1] = 0.

    ### Third row
    A[2, 0] = (-2 * h / Stressobject.tsp) / (2 * A2)
    A[2, 1] = (2 * Stressobject.lsk / Stressobject.tsk + 2 * h / Stressobject.tsp) / (2 * A2)
    A[2, 2] = -1
    b[2] = 0.

    solution = np.linalg.solve(A, b)
    Stressobject.Tq1f = solution[0]
    Stressobject.Tq2f = -solution[0] + solution[1]
    Stressobject.Tq3f = solution[1]
    Stressobject.Tq4f = solution[1]
    Stressobject.Tq5f = -solution[0] + solution[1]
    Stressobject.Tq6f = solution[0]

    ### Auxiliary functions
    Stressobject.compute_stressdistributions(Sy, Sz, My, Mz, T)

    Stressobject.plot_shearflowdistributions()
    Stressobject.plot_directstressdistributions()
    Stressobject.plot_vonmisesstressdistributions()

