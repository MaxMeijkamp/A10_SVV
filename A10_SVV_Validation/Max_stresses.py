#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 15:23:36 2020

@author: mustafawahid
"""
import numpy as np
from validation import *

def node_x(node):
    # returns x-coordinate OF GIVEN NODE
    nodelist = np.genfromtxt("B737.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")
    nodelist = nodelist.astype(np.float)
    for i in nodelist:
        if int(i[0]) == node:
            return i[1]


def nodes(element):
    # enter element number and gives corresponding nodes
    elementlist = np.genfromtxt("B737.inp", dtype=str, skip_header=6598, skip_footer=(14594 - 13233), delimiter=",")
    elementlist = elementlist.astype(np.float)
    for i in elementlist:
        if int(i[0]) == element:
            return i[1], i[2], i[3], i[4]


def elements(node):
    # enter node number and gives corresponding elements
    elementlist = np.genfromtxt("B737.inp", dtype=str, skip_header=6598, skip_footer=(14594 - 13233), delimiter=",")
    elementlist = elementlist.astype(np.float)
    elements = []
    for i in elementlist:
        if int(i[1]) == node:
            elements.append(i[0])
        if int(i[2]) == node:
            elements.append(i[0])
        if int(i[3]) == node:
            elements.append(i[0])
        if int(i[4]) == node:
            elements.append(i[0])
    return elements


def integrationpoint(element):
    elementlist = np.genfromtxt("B737.inp", dtype=str, skip_header=6598, skip_footer=(14594 - 13233), delimiter=",")
    elementlist = elementlist.astype(np.float)
    allnodes = nodes(element)
    x = 0
    for i in allnodes:
        x = x + node_x(i)
    return x / len(allnodes)


# fix this part
def crossection(node):
    # finds all the nodes in the crossection of given node
    elementlist = np.genfromtxt("B737.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")
    elementlist = elementlist.astype(np.float)
    crossection_list = []

    for i in elementlist:
        if int(i[0]) == node:  # find line with corresponding element
            for n in elementlist:
                if int(n[1]) == int(i[1]):  # finds all the elements with the same x coordinate as the input element
                    crossection_list.append(int(n[0]))

    return np.array(crossection_list)


def crossection_element(element):
    a1 = crossection(nodes(element)[0])
    a2 = crossection(nodes(element)[2])
    a = list(np.hstack((a1, a2)))

    lst = []  # list of all the nodes of the crossection of nodes on the left and right of the element
    for i in a:
        lst = lst + elements(i)
    crossectionlist = []
    for i in lst:
        if lst.count(i) == 4:
            if crossectionlist.count(i) == 0:
                crossectionlist.append(i)
        if lst.count(i) == 6:
            if crossectionlist.count(i) == 0:
                crossectionlist.append(i)

    return np.sort(crossectionlist)


def max_von_mises(case):
    # finds the maximum von mises stress at given loading case for the entire aileron
    stresslist = []
    for i in get_dat(case, 'stresses'):
        stresslist.append((i[1], i[0]))
    return max(stresslist)


def max_s12(case):
    # finds maximum s12 at given loading case
    stresslist = []
    for i in get_dat(case, 'stresses'):
        stresslist.append((i[2], i[0]))
    return max(stresslist)


def von_mises(case, element):
    # finds von mises stress at given loading case and element
    for i in get_dat(case, 'stresses'):
        if int(i[0]) == element:
            return i[1]


# enter element number and will give the corresponding S12
def s12(case, element):
    # finds von mises stress at given loading case and element
    for i in get_dat(case, 'stresses'):
        if int(i[0]) == element:
            return i[2]


def crossection_von_mises(case):
    a = max_von_mises(case)[1]
    b = crossection_element(a)
    crossection_von_mises_list = []

    for i in b:
        crossection_von_mises_list.append(von_mises(case, i))

    return crossection_von_mises_list


def crossection_s12(case):
    a = max_s12(case)[1]
    b = crossection_element(a)
    crossection_s12_list = []

    for i in b:
        crossection_s12_list.append(s12(case, i))

    return crossection_s12_list


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
c = get_dat('bending', 'disp')
c = c.astype(np.float)
hingelinedata_bend = []
x_bend = []
y_bend = []
z_bend = []
mag_bend = []

# Deflection data - loading case: jam_bent
d = get_dat('Jam_Bent', 'disp')
d = d.astype(np.float)
hingelinedata_jam_bent = []
x_jam_bent = []
y_jam_bent = []
z_jam_bent = []
mag_jam_bent = []

# Deflection data - loading case: jam_bent
e = get_dat('Jam_Straight', 'disp')
e = e.astype(np.float)
hingelinedata_jam_straight = []
x_jam_straight = []
y_jam_straight = []
z_jam_straight = []
mag_jam_straight = []

# Find hingeline data for loading cases
s = 0
for i in hingelinenodes:
    i = int(i) - 1
    hingelinedata_bend.append(c[i])
    mag_bend.append(c[i][1])
    x_bend.append(c[i][2] + hingelinenodes_x_coordinate[s])
    y_bend.append(c[i][3])
    z_bend.append(c[i][4])

    hingelinedata_jam_bent.append(d[i])
    mag_jam_bent.append(d[i][1])
    x_jam_bent.append(d[i][2] + hingelinenodes_x_coordinate[s])
    y_jam_bent.append(d[i][3])
    z_jam_bent.append(d[i][4])

    hingelinedata_jam_straight.append(e[i])
    mag_jam_straight.append(e[i][1])
    x_jam_straight.append(e[i][2] + hingelinenodes_x_coordinate[s])
    y_jam_straight.append(e[i][3])
    z_jam_straight.append(e[i][4])

    s = s + 1

# plt.title('Vertical Hingeline deflection x-y plane')
# plt.xlabel('Location x-axis[mm]')
# plt.ylabel('Vertical deflection [mm]')
# plt.plot(y_bend, z_bend, 'bo') #bluedots
# plt.plot(y_jam_bent, z_jam_bent, 'ro') #reddots
# plt.plot(y_jam_straight, z_jam_straight, 'go') #greendots
# plt.show()























