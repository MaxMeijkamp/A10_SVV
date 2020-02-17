import numpy as np
from mpl_toolkits import mplot3d  # 3d plotting
import matplotlib.pyplot as plt

# input file
a = np.genfromtxt("B737INP.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")
a = a.astype(np.float)
xin = []
yin = []
zin = []

# output file
b = np.genfromtxt("B737INP.inp", dtype=str, skip_header=14146, skip_footer=(14594 - 14178), delimiter=",", comments="*")
b = b.astype(np.float)

f = open("B737INP.inp", "r")
lines = f.readlines()
f.close()

nodes = list(lines[9:6597])  # node number with x y z coord
elem = np.array(lines[6599:13232])  # elem number with 4 corresponding nodes

for i in a:
    xin.append(i[1])
    yin.append(i[2])
    zin.append(i[3])

fig = plt.figure()
ax = plt.axes(projection='3d')

# Data for elem on aileron
ax.scatter3D(xin, yin, zin)