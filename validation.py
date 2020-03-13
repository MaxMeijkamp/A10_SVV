import numpy as np
from operator import itemgetter
from mpl_toolkits import mplot3d  # 3d plotting
import matplotlib.pyplot as plt

# # input file
#a = np.genfromtxt("A10_SVV_DataSets/B737INP.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")
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
#dd
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
#get the x , y , z of the nodes in diff arrays
# xin = np.genfromtxt("A10_SVV_DataSets/B737INP.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598),
#                   delimiter=",", usecols=1).astype(float)
# yin = np.genfromtxt("A10_SVV_DataSets/B737INP.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598),
#                   delimiter=",", usecols=2).astype(float)
# zin = np.genfromtxt("A10_SVV_DataSets/B737INP.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598),
#                   delimiter=",", usecols=3).astype(float)
# LEa = []  # nmbr, x , y , z
# for n in np.genfromtxt("A10_SVV_DataSets/B737INP.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598),
#                       delimiter=",").astype(float):
#     if n[2] == 0.0 and n[3] == 102.5:
#         LEa.append(n)
# LEa = np.array(LEa)
# print(LEa)
# #plotting the 3D aileron
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.scatter3D(xin, yin, zin)
# ax.scatter3D(LEa[:,1],LEa[:,2],LEa[:,3])
# plt.show()
# print(yin.min(),zin.min())

def get_dat(case, param):

    file = "A10_SVV_DataSets/B737RPT.rpt"
    newload = []

    if str(case) == 'bending' :

        if param == 'stresses': # returns VMS and S12 at each element
            loaded1 = np.genfromtxt(file, dtype=str, skip_header=20, max_rows=5778).astype(float)
            loaded2 = np.genfromtxt(file, dtype=str, skip_header=5816, max_rows=856).astype(float)
            loaded = np.concatenate((loaded1,loaded2))

            for elem in loaded:
                newload.append([elem[0],(elem[2]+elem[3])/2,(elem[4]+elem[5])/2])
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
#print(get_dat('bending','disp'))
#print(get_dat('Jam_Bent','disp'))

H = []

for s in np.genfromtxt("A10_SVV_DataSets/B737INP.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598),
                      delimiter=",").astype(float):
    if s[2]==0 and s[3]==0:
        H.append([s[0],s[1]])

H = H=np.array(H)
print(H)
disph = []

for t in get_dat('bending','disp'):
    #print(t)
    for h in H:
        if int(t[0]) == h[0]:
            disph.append([float(h[1]),float(t[3])])
disph = np.sort(disph, axis=0)

print(disph)
for elem in disph:
    print(elem[1])


# def get_twist(case):
#
#     # get INP data
#
#     a = np.genfromtxt("A10_SVV_DataSets/B737INP.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598),
#                       delimiter=",").astype(float)
#
#     # get rows for specific loading
#
#     if case == 'bending':
#         b = np.genfromtxt("A10_SVV_DataSets/B737RPT.rpt", dtype=str, skip_header=20074, max_rows=6588).astype(float)
#     if case == 'Jam_Bent':
#         b = np.genfromtxt("A10_SVV_DataSets/B737RPT.rpt", dtype=str, skip_header=26724, max_rows=6588).astype(float)
#     if case == 'Jam_Straight':
#         b = np.genfromtxt("A10_SVV_DataSets/B737RPT.rpt", dtype=str, skip_header=33374, max_rows=6588).astype(float)
#
#     # classifying the nodes with their x,y,z locations and the elements with their 4 nodes
#
#     nodes = list(a[:6588])  # node number with x y z coord
#     elem = list(a[6599:13232])  # elem number with 4 corresponding nodes
#
#     # find hingeline nodes
#
#     hingelinenodes = [] # x_position , node nmbr
#     hingeline_x = []
#
#     for i in a:
#         if i[2] == 0:
#             if i[3] == 0:
#                 hingelinenodes.append([i[1],i[0]])
#                 hingeline_x.append(i[1])
#
#
#     # find LE
#
#     LE = []  # LE nodes coord. ( x, number )
#     LEnode = []
#
#     for n in nodes:
#
#         if n[2] == 0.0 and round(n[3],2) == 102.5:  # y = 0 and z = 102.5 are the coordinates of the LE
#             LE.append([n[1],n[0]]) #,n[2],n[3]])
#             LEnode.append(n[0])
#
#
#     LE_LO = []  # LE node displacement LOADED CASE (x position, x disp, y disp)
#     i = 0
#     for LEO in LEnode:
#         for b0 in b:
#             if LEO == int(b0[0]):
#                 # LE_LO.append([b0[0],b0[1],b0[2],b0[3],b0[4]])
#                 LE_LO.append([LE[i][0],b0[2], b0[3]])
#                 i += 1
#     LE_LO = np.sort(LE_LO,axis=0)
#
#     # hinge position after load is applied
#
#     # hingenewf = hingelinenodes # x_position, node nmbr , y disp
#     # hingenew = [] # x_position, y disp
#     # f = 0
#     #
#     # for s in b:
#     #     for t in hingelinenodes:
#     #         if s[0] == t[1]:
#     #             hingenew.append([t[0],s[3]])
#     #             hingenewf[f].append(s[3])
#     #             f += 1
#     # hingenewf.sort()
#     # hingenew.sort()
#
#     #Find TE
#     TE = []  # LE nodes coord. ( x, number )
#     TEnodes = []
#     i = 0
#
#     for n in nodes:
#         if n[2] == 0.0 and round(n[3], 1) == -502.5:  # y = 0 and z = -502.5 are the coordinates of the LE
#             TE.append([n[1], n[0]])  # ,n[2],n[3]])
#             TEnodes.append(n[0])
#
#
#     TEdisp = []
#     num = 0
#     for f in TEnodes:
#         for p in b:
#             if f == p[0]:
#                 TEdisp.append([TE[num][0],p[2],p[3]])
#                 num += 1
#     TE = np.sort(TE,axis=0)
#     for elem in TE:
#         print(TE[i][0]-TE[i+1][0])
#         i = i+1
#     # Twist angle as a function of x (span) [node number, twist angle (rad), y disp (mm)]
#
#     print(len(TE), TE, np.sort(LE,axis=0))
#
#     twistcr = []
#     for samp in range(len(LE_LO)) :
#
#         theta = np.arctan((TEdisp[samp][2]-LE_LO[samp][2])/102.5)
#         #print(LE_LO[samp][3],hingenew[samp])
#         twistcr.append(theta)
#
#     # twistcr.sort()
#     hingeline_x.sort()
#
#     # print(hingelinenodes[:,[0]])
#
#     plt.plot(hingeline_x,twistcr)
#     plt.title('Rate of twist for : '+case)
#     plt.xlabel('x - location')
#     plt.ylabel('Angle [rad]')
#
#     plt.show()
#
#
#
#
#     return twistcr
#
# print(get_twist('Jam_Straight'))



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

    # leading edge nodes
    LE = []  # nmbr, x , y , z
    for n in a[:6588]:
        if n[2] == 0.0 and n[3] == 102.5:
            LE.append(n)
    LE = np.array(LE)

    # trailing edge nodes
    TE = []  # nmbr, x, y, z
    for n in a[:6588]:
        if n[3] == -502.5:
            TE.append(n)

    TE = np.array(TE)

    # Leading edge displacement
    LEdisp = []  # nmbr, magn, dx, dy , dz
    for L in LE:
        for b0 in b:
            if L[0] == b0[0]:
                LEdisp.append(b0)
    LEdisp = np.array(LEdisp)

    # Trailing edge displacement
    TEdisp = []  # nmbr, magn, dx, dy ,dz
    for r in TE:
        for b1 in b:
            if r[0] == b1[0]:
                TEdisp.append(b1)
    TEdisp = np.array(TEdisp)

    # Twist
    TEfull = []  # x , dy, Dz
    for elem in TE:
        for elems in TEdisp:
            if elem[0] == elems[0]:
                TEfull.append([elem[1], elems[3], abs(elems[4] - elem[3])])
    LEfull = []  # x, dy, Dz
    for elem in LE:
        for elems in LEdisp:
            if elem[0] == elems[0]:
                LEfull.append([elem[1], elems[3], abs(elems[4] - elem[3])])
    LEfull = np.sort(LEfull, axis=0)
    TEfull = np.sort(TEfull, axis=0)



    thetas = []
    for i in range(len(LE)):
        theta = np.arctan((LEfull[i][1] - TEfull[i][1]) / (LEfull[i][2] + TEfull[i][2]))
        #theta = np.arctan((LEfull[i][1] - TEfull[i][1]) / (102.5+502.5))
        thetas.append(theta)
        print(theta)

    thetas = np.array(thetas)

    return #np.sort(TE[:,1]), thetas

