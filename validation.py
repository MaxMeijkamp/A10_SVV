import numpy as np
from mpl_toolkits import mplot3d  # 3d plotting
import matplotlib.pyplot as plt

# input file
a = np.genfromtxt("A10_SVV_DataSets/B737INP.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")
a = a.astype(np.float)
xin = []
yin = []
zin = []

# output file
b = np.genfromtxt("A10_SVV_DataSets/B737RPT.rpt", dtype=str, skip_header=20074, skip_footer=(59953-26768))
b = b.astype(np.float)


# classifying the nodes with their x,y,z locations and the elements with their 4 nodes
nodes = list(a[9:6597])  # node number with x y z coord
elem = np.array(a[6599:13232])  # elem number with 4 corresponding nodes

for i in a:
    xin.append(i[1])
    yin.append(i[2])
    zin.append(i[3])

# plotting the 3D aileron

fig = plt.figure()
# ax = plt.axes(projection='3d')

# Data for elem on aileron
# ax.scatter3D(xin, yin, zin)
# plt.show()

# ANGLE OF TWIST
LE = []  # LE nodes coord.
for n in nodes:
    if n[2] == float(0) and n[3] == 102.5:  # y = 0 and z = 102.5 are the coordinates of the LE
        LE.append(n[0])

# print(b)

# LE coord CASE 1 : BENDING
LE_Bend = []  # LE node location BENDING CASE (node number, x,y,z location)
for LEO in LE:
    for b0 in b:
        if LEO == int(b0[0]):
            LE_Bend.append(b0)

# Twist angle as a function of x (span)
twist_case1 = []
angle_case1 = []
for samp in LE_Bend:
    theta = np.arctan(samp[3]/samp[4])
    twist_case1.append(samp[0])
    angle_case1.append(theta)

# ax.scatter3D(LE[1], LE[2], LE[3])

plt.plot(twist_case1, angle_case1)
plt.show()
