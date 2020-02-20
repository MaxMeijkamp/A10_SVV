#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 13:56:00 2020

@author: mustafawahid
"""

import numpy as np
import matplotlib.pyplot as plt

# input file
a = np.genfromtxt("B737.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")
a = a.astype(np.float)


hingelinenodes = []
hingelinenodes_x = []
for i in a:
    if i[2] == 0:
        if i[3] == 0:
            hingelinenodes.append(i[0])
            hingelinenodes_x.append(i[1])

# Deflection data - loading case: bending
c = np.genfromtxt("B737.rpt", skip_header=20074, skip_footer=(59956 - 26662 - 109))
c = c.astype(np.float)
hingelinedata_bend = []
x_bend = []
y_bend = []
z_bend = []
mag_bend = []

# Deflection data - loading case: jam_bent
d = np.genfromtxt("B737.rpt", skip_header=26724, skip_footer=(59956 - 33313 - 89))
d = d.astype(np.float)
hingelinedata_jam_bent = []
x_jam_bent = []
y_jam_bent = []
z_jam_bent = []
mag_jam_bent = []

# Deflection data - loading case: jam_bent
e = np.genfromtxt("B737.rpt", skip_header=33374, skip_footer=(59956 - 39963 - 70))
e = e.astype(np.float)
hingelinedata_jam_straight = []
x_jam_straight = []
y_jam_straight = []
z_jam_straight = []
mag_jam_straight = []

# Find hingeline data for loading cases
for i in hingelinenodes:
    i = int(i) - 1

    hingelinedata_bend.append(c[i])
    mag_bend.append(c[i][1])
    x_bend.append(c[i][2])
    y_bend.append(c[i][3])
    z_bend.append(c[i][4])

    hingelinedata_jam_bent.append(d[i])
    mag_jam_bent.append(d[i][1])
    x_jam_bent.append(d[i][2])
    y_jam_bent.append(d[i][3])
    z_jam_bent.append(d[i][4])

    hingelinedata_jam_straight.append(e[i])
    mag_jam_straight.append(e[i][1])
    x_jam_straight.append(e[i][2])
    y_jam_straight.append(e[i][3])
    z_jam_straight.append(e[i][4])

plt.title('Vertical Hingeline deflection x-y plane')
plt.xlabel('Location x-axis[mm]')
plt.ylabel('Vertical deflection [mm]')
plt.plot(hingelinenodes_x, y_bend, 'o')
plt.plot(hingelinenodes_x, y_jam_bent, 'o')
plt.plot(hingelinenodes_x, y_jam_straight, 'o')
plt.show()