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
    for i in elementlist:
        if int(i[0]) == element:
            integrationpoint_x = (node_x(int(max(allnodes))) + node_x(int(min(allnodes)))) * 0.5
            return integrationpoint_x


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


# nodeslist = np.genfromtxt("B737.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")


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


print (crossection_element(max_von_mises('bending')[1]))
print (crossection_von_mises('bending'))

# element=5190
# nodes_list=nodes(element)
# crossection1 = crossection(min(nodes(element)))
# crossection2 = crossection(max(nodes(element)))
# crossectionlist_nodes=crossection1+crossection2
# crossectionlist=[]
# for i in range (6000):
#    nodes_list=nodes(i)
#    print(nodes_list)


# def crossection_element(element):
#    crossection1 = crossection(min(nodes(element)))
#    crossection2 = crossection(max(nodes(element)))
#    crossectionlist_nodes=crossection1+crossection2
#    crossectionlist=[]
#    for i  in  range(1,max(crossectionlist_nodes)+1):
#        if crossectionlist_nodes.count(nodes(i)[0])+crossectionlist_nodes.count(nodes(i)[1])+crossectionlist_nodes.count(nodes(i)[2])+crossectionlist_nodes.count(nodes(i)[3])=4:
#            crossectionlist.append(i)
#    return crossectionlist
#
# print(crossection_element(5190))

#    a = nodes(i)
#    if crossectionlist_nodes.count(a[0])+crossectionlist_nodes.count(a[1])+crossectionlist_nodes.count(a[2])+crossectionlist_nodes.count(a[3])=4:
#        crossectionlist.append(i)
# print(crossectionlist)


# def crossection_elements(element):
#    crossection1 = crossection(min(nodes(element)))
#    crossection2 = crossection(max(nodes(element)))
#    
#    crossectionlist_nodes=crossection1+crossection2
#    crossectionlist = []

# print (nodes(5190))


# import matplotlib.pyplot as plt
#
## input file
# nodelist = np.genfromtxt("B737.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")
# nodelist = nodelist.astype(np.float)
#
# hingelinenodes = []
# hingelinenodes_x_coordinate = []
# for i in nodelist:
#    if i[2] == 0:
#        if i[3] == 0:
#            hingelinenodes.append(i[0])
#            hingelinenodes_x_coordinate.append(i[1])
#
## Deflection data - loading case: bending
# c = get_dat('bending','disp')
# c = c.astype(np.float)
# hingelinedata_bend = []
# x_bend = []
# y_bend = []
# z_bend = []
# mag_bend = []
# print(c)
## Deflection data - loading case: jam_bent
# d = get_dat('Jam_Bent','disp')
# d = d.astype(np.float)
# hingelinedata_jam_bent = []
# x_jam_bent = []
# y_jam_bent = []
# z_jam_bent = []
# mag_jam_bent = []
# print(d)
## Deflection data - loading case: jam_bent
# e = get_dat('Jam_Straight','disp')
# e = e.astype(np.float)
# hingelinedata_jam_straight = []
# x_jam_straight = []
# y_jam_straight = []
# z_jam_straight = []
# mag_jam_straight = []
# print(e)
## Find hingeline data for loading cases
# s=0
# for i in hingelinenodes:
#    i = int(i) - 1
#    hingelinedata_bend.append(c[i])
#    mag_bend.append(c[i][1])
#    x_bend.append(c[i][2]+hingelinenodes_x_coordinate[s])
#    y_bend.append(c[i][3])
#    z_bend.append(c[i][4])
#
#    hingelinedata_jam_bent.append(d[i])
#    mag_jam_bent.append(d[i][1])
#    x_jam_bent.append(d[i][2]+hingelinenodes_x_coordinate[s])
#    y_jam_bent.append(d[i][3])
#    z_jam_bent.append(d[i][4])
#
#    hingelinedata_jam_straight.append(e[i])
#    mag_jam_straight.append(e[i][1])
#    x_jam_straight.append(e[i][2]+hingelinenodes_x_coordinate[s])
#    y_jam_straight.append(e[i][3])
#    z_jam_straight.append(e[i][4])
#
#    s=s+1
#
# plt.title('Vertical Hingeline deflection x-y plane')
# plt.xlabel('Location x-axis[mm]')
# plt.ylabel('Vertical deflection [mm]')
# plt.plot(x_bend, y_bend, 'bo') #bluedots
# plt.plot(x_jam_bent, y_jam_bent, 'ro') #reddots
# plt.plot(x_jam_straight, y_jam_straight, 'go') #greendots
# plt.show()


