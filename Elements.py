import numpy as np
# input file

def elements(n):
    elementlist = np.genfromtxt("B737.inp", dtype=str, skip_header=6598, skip_footer=(14594 - 13233), delimiter=",")
    elementlist = elementlist.astype(np.float)
    for i in elementlist:
        if int(i[0])==n:
            return i[1],i[2],i[3],i[4]

totalstresslist_bending = np.genfromtxt("B737.rpt", dtype=str, skip_header=20, skip_footer=(59956 - 5799), delimiter=",")
totalstresslist_bending = totalstresslist_bending.astype(np.float)

totalstresslist_jam_bent = np.genfromtxt("B737.rpt", dtype=str, skip_header=6705, skip_footer=(59956 - ), delimiter=",")
totalstresslist_jam_bent = totalstresslist_jam_bent.astype(np.float)

totalstresslist_jam_straight = np.genfromtxt("B737.rpt", dtype=str, skip_header=9, skip_footer=(59956 - 6598), delimiter=",")
totalstresslist_jam_straight = totalstresslist_jam_straight.astype(np.float)