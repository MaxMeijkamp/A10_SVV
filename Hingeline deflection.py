#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 13:56:00 2020

@author: mustafawahid
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

# input file
a = np.genfromtxt("B737.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")
a = a.astype(np.float)

# finding hingeline nodes
lst = []
for i in a:
    if i[2] == 0:
        if i[3] == 0:
            lst.append(i[0])

# find corresponding hingeline location - loading case: bending
c = np.genfromtxt("B737.rpt", skip_header=20074, skip_footer=(59956 - 26662 - 109))
c = c.astype(np.float)
print (c[0])
print (c[-1])

# find corresponding hingeline location - loading case: jam_bent
d = np.genfromtxt("B737.rpt", skip_header=26724, skip_footer=(59956 - 33313 - 89))
d = d.astype(np.float)
print (d[0])
print (d[-1])

# find corresponding hingeline location - loading case: jam_bent
e = np.genfromtxt("B737.rpt", skip_header=33374, skip_footer=(59956 - 39963 - 70))
e = e.astype(np.float)
print (e[0])
print (e[-1])

# hingelinedata=[]
# for i in lst:
#    i=int(i)-1
#    hingelinedata.append(c[i])
#    
# x_bend=[]
# y_bend=[]
# z_bend=[]
# mag_bend=[]
# for i in hingelinedata:
#    mag_bend.append(i[0])
#    x_bend.append(i[1])
#    y_bend.append(i[2])
#    z_bend.append(i[3])


