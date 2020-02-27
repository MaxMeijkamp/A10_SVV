#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 15:23:36 2020

@author: mustafawahid
"""
import numpy as np


def get_dat(case, param):
    file = "B737.rpt"
    newload = []

    if str(case) == 'bending':

        if param == 'stresses':  # returns VMS and S12 at each element
            loaded1 = np.genfromtxt(file, dtype=str, skip_header=20, max_rows=5778).astype(float)
            loaded2 = np.genfromtxt(file, dtype=str, skip_header=5816, max_rows=856).astype(float)
            loaded = np.concatenate((loaded1, loaded2))

            for elem in loaded:
                newload.append([elem[0], (elem[2] + elem[3]) / 2, (elem[4] + elem[5]) / 2])
        if param == 'disp':  # returns displacement at each node
            newload = np.genfromtxt(file, dtype=str, skip_header=20074, max_rows=6588)

    if case == 'Jam_Bent':

        if param == 'stresses':  # returns VMS and S12 at each element
            loaded1 = np.genfromtxt(file, dtype=str, skip_header=6705, max_rows=5778)
            loaded2 = np.genfromtxt(file, dtype=str, skip_header=12501, max_rows=856)
            loaded = np.concatenate((loaded1, loaded2)).astype(float)

            for elem in loaded:
                newload.append([elem[0], (elem[2] + elem[3]) / 2, (elem[4] + elem[5]) / 2])

        if param == 'disp':  # returns displacement at each node
            newload = np.genfromtxt(file, dtype=str, skip_header=26724,
                                    max_rows=6588)

    if case == 'Jam_Straight':
        if param == 'stresses':  # returns VMS and S12 at each element
            loaded1 = np.genfromtxt(file, dtype=str, skip_header=13390, max_rows=5778)
            loaded2 = np.genfromtxt(file, dtype=str, skip_header=19186, max_rows=856)
            loaded = np.concatenate((loaded1, loaded2)).astype(float)

            for elem in loaded:
                newload.append([elem[0], (elem[2] + elem[3]) / 2, (elem[4] + elem[5]) / 2])

        if param == 'disp':  # returns displacement at each node
            newload = np.genfromtxt(file, dtype=str, skip_header=33374, max_rows=6588)

    return (newload)


def elements(n):
    elementlist = np.genfromtxt("B737.inp", dtype=str, skip_header=6598, skip_footer=(14594 - 13233), delimiter=",")
    elementlist = elementlist.astype(np.float)
    for i in elementlist:
        if int(i[0]) == n:
            return i[1], i[2], i[3], i[4]


def nodes(n):
    elementlist = np.genfromtxt("B737.inp", dtype=str, skip_header=6598, skip_footer=(14594 - 13233), delimiter=",")
    elementlist = elementlist.astype(np.float)
    elements = []
    for i in elementlist:
        if int(i[1]) == n:
            elements.append(i[0])
        if int(i[2]) == n:
            elements.append(i[0])
        if int(i[3]) == n:
            elements.append(i[0])
        if int(i[4]) == n:
            elements.append(i[0])
    return elements


# fix this part
def crossection(element):
    # finds all the elements in the crossection of given element
    elementlist = np.genfromtxt("B737.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")
    elementlist = elementlist.astype(np.float)
    crossection_list = []
    for i in elementlist:
        if int(i[0]) == element:  # find line with corresponding element
            for n in elementlist:
                if int(n[1]) == int(i[1]):  # finds all the elements with the same x coordinate as the input element
                    crossection_list.append(int(n[0]))
            return crossection_list


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
    b = crossection(a)
    crossection_von_mises_list = []

    for i in b:
        crossection_von_mises_list.append(von_mises(case, i))

    return crossection_von_mises_list


def crossection_s12(case):
    a = max_s12(case)[1]
    b = crossection(a)
    crossection_s12_list = []

    for i in b:
        crossection_s12_list.append(s12(case, i))

    return crossection_s12_list


print(von_mises('bending', 5190))
print (crossection(5190))

# import matplotlib.pyplot as plt
#
## input file
# a = np.genfromtxt("B737.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")
# a = a.astype(np.float)
#
#
# hingelinenodes = []
# hingelinenodes_x = []
# for i in a:
#    if i[2] == 0:
#        if i[3] == 0:
#            hingelinenodes.append(i[0])
#            hingelinenodes_x.append(i[1])
#
# print(hingelinenodes)
## Deflection data - loading case: bending
# c = np.genfromtxt("B737.rpt", skip_header=20074, skip_footer=(59956 - 26662 - 109))
# c = c.astype(np.float)
# hingelinedata_bend = []
# x_bend = []
# y_bend = []
# z_bend = []
# mag_bend = []
#
## Deflection data - loading case: jam_bent
# d = np.genfromtxt("B737.rpt", skip_header=26724, skip_footer=(59956 - 33313 - 89))
# d = d.astype(np.float)
# hingelinedata_jam_bent = []
# x_jam_bent = []
# y_jam_bent = []
# z_jam_bent = []
# mag_jam_bent = []
#
## Deflection data - loading case: jam_bent
# e = np.genfromtxt("B737.rpt", skip_header=33374, skip_footer=(59956 - 39963 - 70))
# e = e.astype(np.float)
# hingelinedata_jam_straight = []
# x_jam_straight = []
# y_jam_straight = []
# z_jam_straight = []
# mag_jam_straight = []
#
## Find hingeline data for loading cases
# s=0
# for i in hingelinenodes:
#    i = int(i) - 1
#    print(hingelinenodes_x[s])
#    hingelinedata_bend.append(c[i])
#    mag_bend.append(c[i][1])
#    x_bend.append(c[i][2]+hingelinenodes_x[s])
#    y_bend.append(c[i][3])
#    z_bend.append(c[i][4])
#
#    hingelinedata_jam_bent.append(d[i])
#    mag_jam_bent.append(d[i][1])
#    x_jam_bent.append(d[i][2]+hingelinenodes_x[s])
#    y_jam_bent.append(d[i][3])
#    z_jam_bent.append(d[i][4])
#
#    hingelinedata_jam_straight.append(e[i])
#    mag_jam_straight.append(e[i][1])
#    x_jam_straight.append(e[i][2]+hingelinenodes_x[s])
#    y_jam_straight.append(e[i][3])
#    z_jam_straight.append(e[i][4])
#
#    s=s+1
#
## plt.title('Vertical Hingeline deflection x-y plane')
## plt.xlabel('Location x-axis[mm]')
## plt.ylabel('Vertical deflection [mm]')
# plt.plot(x_bend, z_bend, 'bo') #bluedots
## plt.plot(x_jam_bent, y_jam_bent, 'ro') #reddots
# plt.plot(x_jam_straight, y_jam_straight, 'go') #greendots
# plt.show()
