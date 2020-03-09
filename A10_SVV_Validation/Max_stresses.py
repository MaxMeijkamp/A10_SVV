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


def node_x(node):
    # returns x-coordinate OF GIVEN NODE
    nodelist = np.genfromtxt("B737.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")
    nodelist = nodelist.astype(np.float)
    for i in nodelist:
        if int(i[0]) == node:
            return i[1]


def node_y(node):
    # returns y-coordinate OF GIVEN NODE
    nodelist = np.genfromtxt("B737.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")
    nodelist = nodelist.astype(np.float)
    for i in nodelist:
        if int(i[0]) == node:
            return i[2]


def node_z(node):
    # returns z-coordinate OF GIVEN NODE
    nodelist = np.genfromtxt("B737.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")
    nodelist = nodelist.astype(np.float)
    for i in nodelist:
        if int(i[0]) == node:
            return i[3]


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
    y = 0
    z = 0
    for i in allnodes:
        x = x + node_x(i)
        y = y + node_y(i)
        z = z + node_z(i)
    return x / len(allnodes), y / len(allnodes), z / len(allnodes)


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

    return crossectionlist


def crossection_coordinates(element):
    a = crossection_element(element)
    lst = []
    for i in a:
        lst.append(integrationpoint(i))
    return lst


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

print (max_von_mises('bending'))
print (max_von_mises('Jam_Bent'))
print (max_von_mises('Jam_Straight'))
print (max_s12('bending'))
print (max_s12('Jam_Bent'))
print (max_s12('Jam_Straight'))





(0.3617365, 5190.0) #von mises bending
(0.4141685, 2391.0) #von mises Jam Bent
(0.178384, 5369.0) #von mises Jam Straight
(0.07755395000000001, 384.0) #s12 bending
(0.1326395, 2390.0) #s12 Jam Bent
(0.0831562, 5181.0) #s12 Jam Straight
























