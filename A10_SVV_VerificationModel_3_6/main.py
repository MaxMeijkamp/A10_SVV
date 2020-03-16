### Importing required packages
import numpy as np
import Energy
import Stiffness
import Stress
import math as m

######################## Part I - parameters as in assignment #######################################
aircraft = "B737" # Write either A320, F100, CRJ700 or Do228 (bear in mind capitals); this is used for aerodynamic loading
# Ca = 0.484  # m
# la = 1.691  # m
# x1 = 0.149  # m
# x2 = 0.554  # m
# x3 = 1.541  # m
# xa = 0.272  # m
# ha = 0.173  # m
# tsk = 0.0011  # m
# tsp = 0.0025  # m
# tst = 0.0012  # m
# hst = 0.014   # m
# wst = 0.018   # m
# nst = 13  # -
# d1 = 0.00681  # m
# d3 = 0.02030  # m
# theta = m.radians(26)  # rad
# P = 37.9*1000  # N
Ca = 0.605  # m
la = 2.661  # m
x1 = 0.172  # m
x2 = 1.211  # m
x3 = 2.591  # m
xa = 0.35  # m
ha = 0.205  # m
tsk = 0.0011  # m
tsp = 0.0028  # m
tst = 0.0012  # m
hst = 0.016   # m
wst = 0.019   # m
nst = 15  # -
d1 = 0.01154  # m
d3 = 0.01840  # m
theta = m.radians(28)  # rad
P = 97.4*1000  # N


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
totarea = crosssection.totarea            # total cross-section area
yc = crosssection.yc                 # y-coordinate of the centroid
zc = crosssection.zc                 # z-coordinate of the centroid
Iyy = crosssection.Iyy                # moment of inertia about y-axis
Izz = crosssection.Izz                # moment of inertia about z-axis

######################## Part III - Torsional stiffness calculations #######################################
### Primary functions
""" If you desire, you may disable this line, and manually overwrite the values listed between lines 60-62"""
crosssection.compute_shearcenter()   # Run the calculations
crosssection.compute_torsionalstiffness()   # Run the calculations

### Access to important results
"""" If you desire, you can manually overwrite these values. """
ysc = crosssection.ysc                 # y-coordinate of the centroid
zsc = crosssection.zsc                 # z-coordinate of the centroid
_ = crosssection.J                   # torsional constant

######################## Part IV - Deflection calculations #######################################
### Definition of additional parameters
N = 20     # Number of basis functions to use in Rayleigh-Ritz method (total number of coefficients is 3*N)
E = 72.9*10**9       # E-modulus (Pa)
G = 27.1*10**9       # G-modulus (Pa)

### Create the aileron object
""" Merges the cross-sectional properties with the spanwise properties (length and material properties)"""
aileron = Energy.Beam(la,crosssection,N,E,G)

"""Define your boundary conditions; see manual for explanations. The shown boundary conditions are the boundary
conditions for the aileron as described in the assignment."""
aileron.addbcss(x1,0.,-ha/2.,-theta,d1)
# aileron.addbcss(x1,0.,-ha/2.,-theta,0)
aileron.addbcss(x1,0.,-ha/2.,m.pi/2-theta,0)
aileron.addbcss(x2,0.,-ha/2.,0,0)
aileron.addbcss(x2,0.,-ha/2.,m.pi/2,0.)
aileron.addbcss(x3,0.,-ha/2.,-theta,d3)
# aileron.addbcss(x3,0.,-ha/2.,-theta,0)
aileron.addbcss(x3,0.,-ha/2.,m.pi/2-theta,0)
aileron.addbcss(x2-xa/2.,ha/2.,0,m.pi/2.-theta,0)

""""Define your applied loading; see manual for explanations."""
# aileron.addfpl(x2+xa/2.,ha/2.,0,m.pi/2.-theta,-P)

### Primary functions
""" The following line computes the deflections. If you do not want to include the aerodynamic loading, simply write
aileron.compute_deflections(). Note that the aerodynamic loading significantly slows down the program.
If you do want to include the aerodynamic loading, let the variable aircraft (see part I) equal "A320", "F100", "CRJ700",
"Do228", whichever you want to include. 
Note that the name should be spelled exactly as listed above. Note that if the aircraft you write is inconsistent with the
geometry you define at the beginning of this file, the program will not return an error, but will simply produce bogus
results."""
# aileron.compute_deflections("B737") ### Switch aerodynamic loading to the aircraft that is being considered
aileron.compute_deflections()

### Auxiliary functions
"""" A number of auxiliary functions and results are given to you. """

## Simplistic plotting procedures for a first check
#aileron.plotv()             # Plot the deflections in y-direction, its derivative, the bending moment about the z-axis, and the shear force in y.
#aileron.plotw()             # Plot the deflections in z-direction, its derivative, the bending moment about the y-axis, and the shear force in z.
#aileron.plotphi()           # Plot the twist distribution, the torque and the distributed torque.

## For custom post-processing of the solution
x = np.linspace(0,la,num = int(la*1000+1))  # Subsequent functions accept numpy-arrays
# Compute the deflections
v, w, phi = aileron.eval(x)       # Compute the three deflections
v1, w1, phi1= aileron.fdeval(x)     # Compute their their first order derivative
v2, w2, phi2 = aileron.sdeval(x)     # Compute their their second order derivative
v3, w3, phi3 = aileron.tdeval(x)     # Compute their their third order derivative
# Compute the loading
Sy = aileron.Sy(x)               # Compute the shear force in y
Sz = aileron.Sz(x)               # Compute the shear force in z
My = aileron.My(x)               # Compute the moment around the y-axis
Mz = aileron.Mz(x)               # Compute the moment around the z-axis
T = aileron.T(x)                # Compute the torque
tau = aileron.tau(x)              # Compute the distributed torque
loadcase = "_1"
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\v"+loadcase+".csv", v, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\v1"+loadcase+".csv", v1, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\v2"+loadcase+".csv", v2, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\v3"+loadcase+".csv", v3, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\w"+loadcase+".csv", w, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\w1"+loadcase+".csv", w1, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\w2"+loadcase+".csv", w2, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\w3"+loadcase+".csv", w3, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\phi"+loadcase+".csv", phi, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\phi1"+loadcase+".csv", phi1, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\phi2"+loadcase+".csv", phi2, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\phi3"+loadcase+".csv", phi3, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\Sy"+loadcase+".csv", Sy, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\Sz"+loadcase+".csv", Sz, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\My"+loadcase+".csv", My, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\Mz"+loadcase+".csv", Mz, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\T"+loadcase+".csv", T, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\tau"+loadcase+".csv", tau, delimiter=",")
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
x = 1.002
Sy = aileron.Sy(x)
Sz = aileron.Sz(x)
My = aileron.My(x)
Mz = aileron.Mz(x)
T = aileron.T(x)

### Primary functions
""""The following line should never be disabled, as its results are used in the auxiliary functions"""
Stressobject.compute_unitstressdistributions()

### Auxiliary functions
Stressobject.compute_stressdistributions(Sy,Sz,My,Mz,T)

### Some plotting functions
#Stressobject.plot_shearflowdistributions()
#Stressobject.plot_directstressdistributions()
#Stressobject.plot_vonmisesstressdistributions()

### Access to important results
theta = np.linspace(0,m.pi/2,num = 100)
a = Stressobject.q1f(theta)             # Compute the shear flow distribution in region 1
b = Stressobject.sigma1f(theta)         # Compute the direct stress distribution in region 1
c1 = Stressobject.vm1(theta)             # Compute the Von Mises stress distribution in region 1
d1, e1 = Stressobject.coord1(theta)       # Compute the z,y-coordinates for region 1

y = np.linspace(ha/2.,ha,num = 100)
_ = Stressobject.q2f(y)             # Compute the shear flow distribution in region 3
_ = Stressobject.sigma2f(y)         # Compute the direct stress distribution in region 3
c2 = Stressobject.vm2(y)             # Compute the Von Mises stress distribution in region 3
d2, e2 = Stressobject.coord2(y)       # Compute the z,y-coordinates for region 3

s = np.linspace(0,m.sqrt((Ca-ha/2.)**2+(ha/2.)**2),num = 100)
_ = Stressobject.q3f(s)             # Compute the shear flow distribution in region 4
_ = Stressobject.sigma3f(s)         # Compute the direct stress distribution in region 4
c3 = Stressobject.vm3(s)             # Compute the Von Mises stress distribution in region 4
d3, e3 = Stressobject.coord3(s)       # Compute the z,y-coordinates for region 4

s = np.linspace(0,m.sqrt((Ca-ha/2.)**2+(ha/2.)**2),num = 100)
_ = Stressobject.q4f(s)             # Compute the shear flow distribution in region 4
_ = Stressobject.sigma4f(s)         # Compute the direct stress distribution in region 4
c4 = Stressobject.vm4(s)             # Compute the Von Mises stress distribution in region 4
d4, e4 = Stressobject.coord4(s)       # Compute the z,y-coordinates for region 4

y = np.linspace(0.,ha/2.,num = 100)
_ = Stressobject.q5f(y)             # Compute the shear flow distribution in region 5
_ = Stressobject.sigma5f(y)         # Compute the direct stress distribution in region 5
c5 = Stressobject.vm5(y)             # Compute the Von Mises stress distribution in region 5
d5, e5 = Stressobject.coord5(y)       # Compute the z,y-coordinates for region 5

theta = np.linspace(-m.pi/2,0,num = 100)
_ = Stressobject.q6f(theta)             # Compute the shear flow distribution in region 6
_ = Stressobject.sigma6f(theta)         # Compute the direct stress distribution in region 6
c6 = Stressobject.vm6(theta)             # Compute the Von Mises stress distribution in region 6
d6, e6 = Stressobject.coord6(theta)       # Compute the z,y-coordinates for region 6

np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\c1"+loadcase+".csv", c1, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\d1"+loadcase+".csv", d1, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\e1"+loadcase+".csv", e1, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\c2"+loadcase+".csv", c2, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\d2"+loadcase+".csv", d2, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\e2"+loadcase+".csv", e2, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\c3"+loadcase+".csv", c3, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\d3"+loadcase+".csv", d3, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\e3"+loadcase+".csv", e3, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\c4"+loadcase+".csv", c4, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\d4"+loadcase+".csv", d4, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\e4"+loadcase+".csv", e4, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\c5"+loadcase+".csv", c5, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\d5"+loadcase+".csv", d5, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\e5"+loadcase+".csv", e5, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\c6"+loadcase+".csv", c6, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\d6"+loadcase+".csv", d6, delimiter=",")
np.savetxt("D:\\Documents\\SVV\\A10\\datafiles2\\e6"+loadcase+".csv", e6, delimiter=",")

# cmax = []
# for x in np.linspace(0,la,int(la*1000+1)):
#     Sy = aileron.Sy(x)
#     Sz = aileron.Sz(x)
#     My = aileron.My(x)
#     Mz = aileron.Mz(x)
#     T = aileron.T(x)
#     Stressobject.compute_unitstressdistributions()
#     Stressobject.compute_stressdistributions(Sy,Sz,My,Mz,T)
#     c1 = Stressobject.vm1(theta)
#     c2 = Stressobject.vm2(theta)
#     c3 = Stressobject.vm3(theta)
#     c4 = Stressobject.vm4(theta)
#     c5 = Stressobject.vm5(theta)
#     c6 = Stressobject.vm6(theta)
#     cs = np.vstack((c1, c2, c3, c4, c5, c6))
#     cmax.append(np.amax(cs))
# print(cmax)
# m = max(cmax)
# print(m)
# print([i for i, j in enumerate(cmax) if j == m])
#5656754447.685532, 1002
#18873195687.441074, 0
#15389112052.809385, 0




