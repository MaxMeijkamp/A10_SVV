import numpy as np

def elements(n):
    elementlist = np.genfromtxt("B737.inp", dtype=str, skip_header=6598, skip_footer=(14594 - 13233), delimiter=",")
    elementlist = elementlist.astype(np.float)
    for i in elementlist:
        if int(i[0])==n:
            return i[1],i[2],i[3],i[4]

#total stresses for every loadcases
totalstresslist_bending = np.genfromtxt("B737.rpt", dtype=str, skip_header=20, skip_footer=(59956 - 5799-165))
totalstresslist_jam_bent = np.genfromtxt("B737.rpt", dtype=str, skip_header=6705, skip_footer=(59956 - 12484-146))
totalstresslist_jam_straight = np.genfromtxt("B737.rpt", dtype=str, skip_header=13390, skip_footer=(59956 - 19169-127))


#von misses and S12 stresses per loadcase and corresponding elements

von_misses_bending = []
s12_bending = []
von_misses_jam_bent = []
s12_jam_bent = []
von_misses_jam_straight = []
s12_jam_straight = []

for i in totalstresslist_bending:

    von_misses_bending.append(((float(i[2]) + float(i[3]))*0.5,int(i[0])))
    s12_bending.append(((float(i[4]) + float(i[5]))*0.5,int(i[0])))

for i in totalstresslist_jam_bent:

    von_misses_jam_bent.append(((float(i[2]) + float(i[3]))*0.5,int(i[0])))
    s12_jam_bent.append(((float(i[4]) + float(i[5]))*0.5,int(i[0])))

for i in totalstresslist_jam_straight:

    von_misses_jam_straight.append(((float(i[2]) + float(i[3]))*0.5,int(i[0])))
    s12_jam_straight.append(((float(i[4]) + float(i[5]))*0.5,int(i[0])))

#enter element number and will give the corresponding von mises stress
def von_mises(element):
    von_mises_list=[]
    for i in totalstresslist_bending:
        if int(i[0]) == int(element):
            von_mises_list.append((float(i[2]) + float(i[3]))*0.5)
    for i in totalstresslist_jam_bent:
        if int(i[0]) == element:
            von_mises_list.append((float(i[2]) + float(i[3])) * 0.5)
    for i in totalstresslist_jam_straight:
        if int(i[0]) == element:
            von_mises_list.append((float(i[2]) + float(i[3])) * 0.5)
    return von_mises_list

#enter element number and will give the corresponding S12
def s12(element):
    S12_list = []
    for i in totalstresslist_bending:
        if int(i[0]) == int(element):
            S12_list.append((float(i[4]) + float(i[5])) * 0.5)
    for i in totalstresslist_jam_bent:
        if int(i[0]) == element:
            S12_list.append((float(i[4]) + float(i[5])) * 0.5)
    for i in totalstresslist_jam_straight:
        if int(i[0]) == element:
            S12_list.append((float(i[4]) + float(i[5])) * 0.5)
    return S12_list

#finds the elements in the crossection for entered element
def crossection(element):
    elementlist = np.genfromtxt("B737.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")
    elementlist = elementlist.astype(np.float)
    crossection_list=[]
    for i in elementlist:
        if int(i[0])==element:
            for n in elementlist:
                if int(n[1])==int(i[1]):
                    crossection_list.append(n[0])

    return crossection_list
vonmiseslistbending=[]
a = max(von_misses_bending)[1]
b = crossection(a)
for i in b:
    vonmiseslistbending.append(von_mises(i))
print (max(vonmiseslistbending))

