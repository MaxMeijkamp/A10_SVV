import numpy as np

def elements(n):
    elementlist = np.genfromtxt("B737.inp", dtype=str, skip_header=6598, skip_footer=(14594 - 13233), delimiter=",")
    elementlist = elementlist.astype(np.float)
    for i in elementlist:
        if int(i[0])==n:
            return i[1],i[2],i[3],i[4]


totalstresslist_bending = np.genfromtxt("B737.rpt", dtype=str, skip_header=20, skip_footer=(59956 - 5799-165))
totalstresslist_jam_bent = np.genfromtxt("B737.rpt", dtype=str, skip_header=6705, skip_footer=(59956 - 12484-146))
totalstresslist_jam_straight = np.genfromtxt("B737.rpt", dtype=str, skip_header=13390, skip_footer=(59956 - 19169-127))

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

print (max(von_misses_bending) ,max(von_misses_jam_bent), max(von_misses_jam_straight))
print (max(s12_bending), max(s12_jam_bent), max(s12_jam_straight))



