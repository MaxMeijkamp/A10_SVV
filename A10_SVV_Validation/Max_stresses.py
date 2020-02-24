import numpy as np

totalstresslist_bending = np.genfromtxt("B737.rpt", dtype=str, skip_header=20, skip_footer=(59956 - 5799-165))
totalstresslist_jam_bent = np.genfromtxt("B737.rpt", dtype=str, skip_header=6705, skip_footer=(59956 - 12484-146))
totalstresslist_jam_straight = np.genfromtxt("B737.rpt", dtype=str, skip_header=13390, skip_footer=(59956 - 19169-127))


von_misses_bending=[]
# print(int(totalstresslist_bending[0][2])+int(totalstresslist_bending[0][2]))
for i in totalstresslist_bending:
    von_misses_bending.append((float(i[2]) + float(i[3]))*0.5)
print (max(von_misses_bending))

von_misses_jam_straight=[]
von_misses_jam_bent=[]
s12_bending=[]
s12_jam_bent=[]
s12_jam_straight=[]

