#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 13:56:00 2020

@author: mustafawahid
"""

import matplotlib.pyplot as plt

# input file
nodelist = np.genfromtxt("B737.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")
nodelist = nodelist.astype(np.float)

hingelinenodes = []
hingelinenodes_x_coordinate = []
for i in nodelist:
    if i[2] == 0:
        if i[3] == 0:
            hingelinenodes.append(i[0])
            hingelinenodes_x_coordinate.append(i[1])

# Deflection data - loading case: bending
c = get_dat('bending','disp')
c = c.astype(np.float)
hingelinedata_bend = []
x_bend = []
y_bend = []
z_bend = []
mag_bend = []

# Deflection data - loading case: jam_bent
d = get_dat('Jam_Bent','disp')
d = d.astype(np.float)
hingelinedata_jam_bent = []
x_jam_bent = []
y_jam_bent = []
z_jam_bent = []
mag_jam_bent = []

# Deflection data - loading case: jam_bent
e = get_dat('Jam_Straight','disp')
e = e.astype(np.float)
hingelinedata_jam_straight = []
x_jam_straight = []
y_jam_straight = []
z_jam_straight = []
mag_jam_straight = []

# Find hingeline data for loading cases
s=0
for i in hingelinenodes:
    i = int(i) - 1
    hingelinedata_bend.append(c[i])
    mag_bend.append(c[i][1])
    x_bend.append(c[i][2]+hingelinenodes_x_coordinate[s])
    y_bend.append(c[i][3])
    z_bend.append(c[i][4])

    hingelinedata_jam_bent.append(d[i])
    mag_jam_bent.append(d[i][1])
    x_jam_bent.append(d[i][2]+hingelinenodes_x_coordinate[s])
    y_jam_bent.append(d[i][3])
    z_jam_bent.append(d[i][4])

    hingelinedata_jam_straight.append(e[i])
    mag_jam_straight.append(e[i][1])
    x_jam_straight.append(e[i][2]+hingelinenodes_x_coordinate[s])
    y_jam_straight.append(e[i][3])
    z_jam_straight.append(e[i][4])

    s=s+1

plt.title('Vertical Hingeline deflection x-y plane')
plt.xlabel('Location x-axis[mm]')
plt.ylabel('Vertical deflection [mm]')
plt.plot(y_bend, z_bend, 'bo') #bluedots
plt.plot(y_jam_bent, z_jam_bent, 'ro') #reddots
plt.plot(y_jam_straight, z_jam_straight, 'go') #greendots
plt.show()



#x1 = np.savetxt("Hingelinedata deflection Bending X location",x_bend,delimiter=",")
#x2 = np.savetxt("Hingelinedata deflection Jam Bent X location",x_jam_bent,delimiter=",")
#x3 = np.savetxt("Hingelinedata deflection Jam straight X location",x_jam_straight,delimiter=",")
#y1 = np.savetxt("Hingelinedata deflection Bending y location",y_bend,delimiter=",")
#y2 = np.savetxt("Hingelinedata deflection Jam Bent y location",y_jam_bent,delimiter=",")
#y3 = np.savetxt("Hingelinedata deflection Jam straight y location",y_jam_straight,delimiter=",")

















