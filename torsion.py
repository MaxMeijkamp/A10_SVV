import numpy as np
from math import *
from InputClasses import *
import matplotlib.pyplot as plt
import os

a1 = Aileron()
##########################################################

d, shear_list, s_56, idx_s, iki = a1.shearcentre()

def num_twist(x, shear_list, s_56, idx_s):
    a = Aileron()

    L = x
    G = 27*10**9
    Am1 = (pi*(a.height/2)**2)/2
    Am2 = ((a.chord-a.height/2)*a.height)/2
    C1 = 1/(2*G*Am1)
    C2 = 1/(2*G*Am2)

    s1 = pi*a.height/2
    s2 = a.height
    s3 = sqrt((a.chord-a.height)**2+(a.height/2)**2)
    s4 = s3
    s5 = s2

    B1 = s1/a.skint + s2/a.spart
    B2 = s2/a.spart
    D1 = (s3+s4)/a.skint + s5/a.spart

    qs02 = 1 / 2*( Am2 + Am1*(C2*D1+B2*C1)/(C1*B1+C2*B2))
    qs01 = (qs02*(C2*D1+B2*C1))/(C1*B1+C2*B2)

    # theta_1 = L * C1 * (qs01*B1 - qs02*B2)
    # theta_2 = L * C2 * (qs02*D1 - qs01*B2)


    #good one


    qs1 = shear_list[3,0:idx_s[1]]
    n1 = len(qs1)
    t1 = a.skint

    qs2 = s_56[3]
    n2 = len(qs2)
    t2 = a.spart

    qs3 = shear_list[3,idx_s[3]:idx_s[4]]
    n3 = len(qs3)

    qs4 = shear_list[3,idx_s[1]:idx_s[2]]
    n4 = len(qs4)

    qs5 = shear_list[3,idx_s[2]:idx_s[3]]
    n5 = len(qs5)

    elem_width = a._circumference / shear_list.size

    X = elem_width


    theta_01 = L*C1*((sum(qs1)-n1*qs01)*n1/t1 - (sum(qs2)+n2*(qs01-qs02))*n2/t2 + (sum(qs3)-n2*qs01)*n3/t1)*X

    theta_02 = L*C2*((sum(qs4)-n4*qs02)*n4/t2 + (sum(qs5)-n5*qs02)*n5/t1 + (sum(qs2)-n2*(qs02-qs01)*n2)/t2)*X


    return theta_01, theta_02,




#list = num_twist()


ang = []
x1 = []

for x in np.linspace(0,a1.span,1000):
    ang.append([num_twist(x, shear_list, s_56, idx_s)[0],num_twist(x, shear_list, s_56, idx_s)[1]])
    x1.append(x)

#print(x1,ang[0:])
plt.plot(x1,[i[0] for i in ang])
plt.plot(x1,[i[1] for i in ang])
plt.show()

if __name__ == "__main__":
    main_path = os.path.dirname(os.path.realpath(__file__))  # Script folder
    data_path = os.path.join(*[main_path, 'datafiles'])  # Data folder


    def dataimport(file):
        return np.genfromtxt(data_path + "\\" + file + ".csv", delimiter=",")

    Sy = dataimport("shearsolved_y")


