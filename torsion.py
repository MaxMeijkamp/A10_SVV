import numpy as np
from math import *
from InputClasses import *

def num_twist(x):
    a = Aileron()
    L = x
    G = 27*10*9
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

    q2 = 1 / 2*( Am2 + Am1*(C2*D1+B2*C1)/(C1*B1+C2*B2))
    q1 = (q2*(C2*D1+B2*C1))/(C1*B1+C2*B2)

    theta_1 = L * C1 * (q1*B1 - q2*B2)
    theta_2 = L * C2 * (q2*D1 - q1*B2)

    return theta_1, theta_2

