import numpy as np
from operator import itemgetter
from mpl_toolkits import mplot3d  # 3d plotting
import matplotlib.pyplot as plt

# # input file
# a = np.genfromtxt("A10_SVV_DataSets/B737INP.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")
# #a = a.astype(np.float)
# xin = []
# yin = []
# zin = []
#
# # output file
# b = np.genfromtxt("A10_SVV_DataSets/B737RPT.rpt", dtype=str, skip_header=20074, skip_footer=(59953-26768))
# b = b.astype(np.float)
#
#
#
# # classifying the nodes with their x,y,z locations and the elements with their 4 nodes
# nodes = list(a[9:6597])  # node number with x y z coord
# elem = np.array(a[6599:13232])  # elem number with 4 corresponding nodes
#
# for i in a:
#     xin.append(i[1])
#     yin.append(i[2])
#     zin.append(i[3])
#
# # plotting the 3D aileron
#
# fig = plt.figure()
# #ax = plt.axes(projection='3d')
#
# # Data for elem on aileron
# # ax.scatter3D(xin, yin, zin)
# # plt.show()
#
# # ANGLE OF TWIST
# LE = [] # LE nodes coord.
# for n in nodes:
#     if n[2] == float(0) and n[3] == 102.5:  # y = 0 and z = 102.5 are the coordinates of the LE
#          LE.append(n[0])
#
# #print(b)
#
# # LE coord CASE 1 : BENDING
# LE_Bend = [] # LE node location BENDING CASE (node number, x,y,z location)
# for LEO in LE:
#     for b0 in b:
#         if LEO == int(b0[0]):
#             LE_Bend.append(b0)
#
# # Twist angle as a function of x (span)
# LE_nodes = []
# angle_case1 = []
# LE_bend1 = []
# for samp in LE_Bend:
#     theta = np.arctan(samp[3]/samp[4])
#     print(samp[4])
#     LE_nodes.append(samp[0])
#     angle_case1.append(theta)
#     LE_bend1.append(samp[3])
#
# # print(LE_Bend)
# # print(len(LE_Bend))
#
#
#
# #ax.scatter3D(LE[1], LE[2], LE[3])
#
# # plt.plot(LE_nodes,angle_case1)
# # plt.plot(LE_nodes, LE_bend1)
# # plt.show()


# plot the aileron
# get the x , y , z of the nodes in diff arrays
#     xin = np.genfromtxt("A10_SVV_DataSets/B737INP.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598),
#                       delimiter=",", usecols=1).astype(float)
#     yin = np.genfromtxt("A10_SVV_DataSets/B737INP.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598),
#                       delimiter=",", usecols=2).astype(float)
#     zin = np.genfromtxt("A10_SVV_DataSets/B737INP.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598),
#                       delimiter=",", usecols=3).astype(float)
#
#     # plotting the 3D aileron
#     # fig = plt.figure()
#     # ax = plt.axes(projection='3d')
#     # ax.scatter3D(xin, yin, zin)
#     # plt.show()
def get_dat(case, param):

    file = "A10_SVV_DataSets/B737RPT.rpt"
    newload = []

    if str(case) == 'bending' :

        if param == 'stresses': # returns VMS and S12 at each element
            loaded1 = np.genfromtxt(file, dtype=str, skip_header=20, max_rows=5778).astype(float)
            loaded2 = np.genfromtxt(file, dtype=str, skip_header=5816, max_rows=856).astype(float)
            loaded = np.concatenate((loaded1,loaded2))

            for elem in loaded:
                newload = newload.append([elem[0],(elem[2]+elem[3])/2,(elem[4]+elem[5])/2])
        if param == 'disp': # returns displacement at each node
            newload = np.genfromtxt(file, dtype=str, skip_header=20074,max_rows = 6588)

    if case == 'Jam_Bent':

        if param == 'stresses':  # returns VMS and S12 at each element
            loaded1 = np.genfromtxt(file, dtype=str, skip_header=6705, max_rows=5778)
            loaded2 = np.genfromtxt(file, dtype=str, skip_header=12501, max_rows=856)
            loaded = np.concatenate((loaded1, loaded2)).astype(float)

            for elem in loaded:
                newload.append([elem[0], (elem[2] + elem[3]) / 2, (elem[4] + elem[5]) / 2])

        if param == 'disp': # returns displacement at each node
            newload = np.genfromtxt(file, dtype=str, skip_header=26724,
                                        max_rows=6588)

    if case == 'Jam_Straight':
        if param == 'stresses':  # returns VMS and S12 at each element
            loaded1 = np.genfromtxt(file, dtype=str, skip_header=13390, max_rows=5778)
            loaded2 = np.genfromtxt(file, dtype=str, skip_header=19186, max_rows=856)
            loaded = np.concatenate((loaded1, loaded2)).astype(float)

            for elem in loaded:
                newload.append([elem[0], (elem[2] + elem[3]) / 2, (elem[4] + elem[5]) / 2])

        if param == 'disp': # returns displacement at each node
            newload = np.genfromtxt(file, dtype=str, skip_header=33374,
                                        max_rows=6588)


    return (newload)



def get_twist(case):

    # get INP data

    a = np.genfromtxt("A10_SVV_DataSets/B737INP.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598),
                      delimiter=",").astype(float)

    # get rows for specific loading

    if case == 'bending':
        b = np.genfromtxt("A10_SVV_DataSets/B737RPT.rpt", dtype=str, skip_header=20074, max_rows=6588).astype(float)
    if case == 'Jam_Bent':
        b = np.genfromtxt("A10_SVV_DataSets/B737RPT.rpt", dtype=str, skip_header=26724, max_rows=6588).astype(float)
    if case == 'Jam_Straight':
        b = np.genfromtxt("A10_SVV_DataSets/B737RPT.rpt", dtype=str, skip_header=33374, max_rows=6588).astype(float)

    # classifying the nodes with their x,y,z locations and the elements with their 4 nodes

    nodes = list(a[:6588])  # node number with x y z coord
    elem = list(a[6599:13232])  # elem number with 4 corresponding nodes

    # find hingeline nodes

    hingelinenodes = []
    hingelinenodes_x = []

    for i in a:
        if i[2] == 0:
            if i[3] == 0:
                hingelinenodes.append(i[0])
                hingelinenodes_x.append(i[1])

    # find LE

    LE = []  # LE nodes coord.


    for n in nodes:
        if n[2] == 0.0 and round(n[3],2) == 102.5:  # y = 0 and z = 102.5 are the coordinates of the LE
            #LE.append([n[0],n[1],n[2],n[3]])
            LE.append(n)

    np.sort(LE, axis=1)


    LE_LO = []  # LE node displacement LOADED CASE (node number, magnitude , x,y,z location)

    for LEO in LE:
        for b0 in b:
            if LEO[0] == int(b0[0]):
                LE_LO.append([b0[0],b0[1],b0[2],b0[3],b0[4]])


    # hinge position after load is applied

    hingenew = [] # node nmbr , y disp

    for s in b:
        for t in hingelinenodes:
            if s[0] == t:
                hingenew.append(s[3])


    # Twist angle as a function of x (span) [node number, twist angle (rad), y disp (mm)]
    twistcr = []
    for samp in range(len(LE_LO)) :

        theta = np.arctan((-hingenew[samp]+LE_LO[samp][3])/102.5)
        #print(LE_LO[samp][3],hingenew[samp])
        twistcr.append(theta)



    twistcr.sort()
    hingelinenodes_x.sort()



    # plt.plot(hingelinenodes_x,twistcr)
    # plt.show()


    return(LE)