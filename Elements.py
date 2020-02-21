import numpy as np
# input file
a = np.genfromtxt("B737.inp", dtype=str, skip_header=6598, skip_footer=(14594 - 13233), delimiter=",")
a = a.astype(np.float)


print (a)
# def elements(n):



