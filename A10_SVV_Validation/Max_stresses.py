import numpy as np

#total stresses for every loadcases and region
#
# totalstresslist_bending_region1 = np.genfromtxt("B737.rpt", dtype=str, skip_header=20, skip_footer=(59956 - 5799-165))
# totalstresslist_bending_region2 = np.genfromtxt("B737.rpt", dtype=str, skip_header=5816, skip_footer=(59956 - 6637 -194))
# totalstresslist_jam_bent_region1 = np.genfromtxt("B737.rpt", dtype=str, skip_header=6705, skip_footer=(59956 - 12484-146))
# totalstresslist_jam_bent_region2 = np.genfromtxt("B737.rpt", dtype=str, skip_header=12501, skip_footer=(59956 - 13358-139))
# totalstresslist_jam_straight_region1 = np.genfromtxt("B737.rpt", dtype=str, skip_header=13390, skip_footer=(59956 - 19169-127))
# totalstresslist_jam_straight_region2 = np.genfromtxt("B737.rpt", dtype=str, skip_header=19186, skip_footer=(59956 - 20043-120))
# totalstresslist_bending=[]




# totalstresslist_bending = np.genfromtxt("B737.rpt", dtype=str, skip_header=20, skip_footer=(59956 - 5799-165))
# totalstresslist_jam_bent = np.genfromtxt("B737.rpt", dtype=str, skip_header=6705, skip_footer=(59956 - 12484-146))
# totalstresslist_jam_straight = np.genfromtxt("B737.rpt", dtype=str, skip_header=13390, skip_footer=(59956 - 19169-127))

#von misses and S12 stresses per loadcase and corresponding elements

# von_misses_bending = []
# s12_bending = []
# von_misses_jam_bent = []
# s12_jam_bent = []
# von_misses_jam_straight = []
# s12_jam_straight = []

# for i in totalstresslist_bending:
#
#     von_misses_bending.append(((float(i[2]) + float(i[3]))*0.5,int(i[0])))
#     s12_bending.append(((float(i[4]) + float(i[5]))*0.5,int(i[0])))
#
# for i in totalstresslist_jam_bent:
#
#     von_misses_jam_bent.append(((float(i[2]) + float(i[3]))*0.5,int(i[0])))
#     s12_jam_bent.append(((float(i[4]) + float(i[5]))*0.5,int(i[0])))
#
# for i in totalstresslist_jam_straight:
#
#     von_misses_jam_straight.append(((float(i[2]) + float(i[3]))*0.5,int(i[0])))
#     s12_jam_straight.append(((float(i[4]) + float(i[5]))*0.5,int(i[0])))


# !/usr/bin/env python3
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


# fix this part
def crossection(element):
    # finds all the nodes in the crossection of given node
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