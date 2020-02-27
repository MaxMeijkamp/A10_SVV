import numpy as np
from math import *
from InputClasses import *


a = Aileron()

##########################################################


def shearcentre(a):
    r = a.radius
    a = a.a
    Izz = a.Izz()
    skint = a.skint

    dist_s1 = dist_s4 = .5 * np.pi * r
    dist_s2 = dist_s3 = a
    dist_s5 = dist_s6 = r

    def func_s(s):
        if s <= dist_s1:
            return (r * np.sin(s / r)) * skint
        elif dist_s1 < s <= dist_s1 + dist_s2:
            return (r - (s - dist_s1) * r / a) * skint
        elif dist_s2 + dist_s1 < s <= dist_s3 + dist_s2 + dist_s1:
            return (0. - (s - dist_s1 - dist_s2) * r / a) * skint
        elif dist_s3 + dist_s2 + dist_s1 < s <= dist_s4 + dist_s3 + dist_s2 + dist_s1:
            return (r * -np.cos((s - .5 * r * np.pi - a * 2) / r)) * skint

    def func_s56(s):  # from s = 0 to s = height / 2
        return s * a.spart

    integrateV = np.vectorize(integrate)
    func_s = np.vectorize(func_s)
    num_steps = 10000

    s_list = np.zeros((7, num_steps))  # array with s value, integral shear, booms
    s_list[0] = np.linspace(0, a._circumference, num_steps)
    idx_s = [
        0, int((np.abs(s_list[0] - dist_s1)).argmin()),
        int(np.abs(s_list[0] - dist_s1 - dist_s2).argmin()),
        int(np.abs(s_list[0] - dist_s1 - dist_s2 - dist_s3).argmin()),
        int(np.abs(s_list[0] - dist_s1 - dist_s2 - dist_s3 - dist_s4).argmin())
    ]

    ss = np.arange(-r, r, a._circumference / num_steps)
    s56_list = np.zeros((7, ss.size))

    s56_list[0] = ss
    idx_s56 = [int(ss.size / 2), int(ss.size - 1)]
    s56_list[1, idx_s56[0]: idx_s56[1]] = - 1. / Izz * (
        integrateV(func_s56, 0, s56_list[0, idx_s56[0]: idx_s56[1]], 100))
    s56_list[1, 0: idx_s56[0]] = - 1. / Izz * (integrateV(func_s56, 0, s56_list[0, 0: idx_s56[0]], 100))

    # Base shear flow of skin:
    s_list[1, 0: idx_s[1]] = -1. / Izz * (integrateV(func_s, 0, s_list[0, 0: idx_s[1]], 100))
    s_list[1, idx_s[1]: idx_s[2]] = -1. / Izz * (
        integrateV(func_s, s_list[0, idx_s[1]], s_list[0, idx_s[1]: idx_s[2]], 100)) + s_list[1, idx_s[1] - 1]
    s_list[1, idx_s[2]: idx_s[3]] = -1. / Izz * (
        integrateV(func_s, s_list[0, idx_s[2]], s_list[0, idx_s[2]: idx_s[3]], 100)) + s_list[1, idx_s[2] - 1]
    s_list[1, idx_s[3]: idx_s[4]] = -1. / Izz * (
        integrateV(func_s, s_list[0, idx_s[3]], s_list[0, idx_s[3]: idx_s[4]], 100)) + s_list[1, idx_s[3] - 1]
    s_list[1, idx_s[4]] = -1. / Izz * (integrateV(func_s, s_list[0, idx_s[3]], s_list[0, idx_s[4]], 100)) + s_list[
        1, idx_s[3] - 1]

    # Base shear flow of booms:
    stepsize = a._circumference / a.stiffn

    for i in range(a.stiffn):
        boomdist = stepsize * i
        idx = (np.abs(s_list[0] - boomdist)).argmin()
        s_list[2, idx: -1] += -1. / Izz * a.stiffener_area * func_s(s_list[0, idx]) / a.skint

    # Combined in row 3:
    s_list[3] = s_list[2] + s_list[1]

    stepcorrection = (num_steps - 1) / num_steps
    # Find correction shear flow:
    self1 = a
    stepwidth = self1._circumference / num_steps
    shearbasesum1 = stepwidth * (
            (np.sum(s_list[3, idx_s[0]: idx_s[1]]) + np.sum(s_list[3, idx_s[3]: idx_s[4]])) / self1.skint - np.sum(
            s56_list[1]) / self1.spart)
    shearbasesum2 = stepwidth * (
            (np.sum(s_list[3, idx_s[1]: idx_s[2]]) + np.sum(s_list[3, idx_s[2]: idx_s[3]])) / self1.skint + np.sum(
            s56_list[1]) / self1.spart)
    twistmatrix = np.array([[2 * dist_s1 / self1.skint + 2 * dist_s5 / self1.spart, -2 * dist_s5 / self1.spart],
                            [-2 * dist_s5 / self1.spart, 2 * dist_s2 / self1.skint + 2 * dist_s5 / self1.spart]])
    shearbasevector = np.array([[-shearbasesum1], [-shearbasesum2]])
    correctionshears = np.linalg.solve(twistmatrix, shearbasevector)
    correctionshear1 = correctionshears[0]
    correctionshear2 = correctionshears[1]

    # Add correction shears and base shear to row 4
    s_list[4, idx_s[0]: idx_s[1]] = s_list[3, idx_s[0]: idx_s[1]] + correctionshear1
    s_list[4, idx_s[3]: idx_s[4]] = s_list[3, idx_s[3]: idx_s[4]] + correctionshear1
    s_list[4, idx_s[1]: idx_s[2]] = s_list[3, idx_s[1]: idx_s[2]] + correctionshear2
    s_list[4, idx_s[2]: idx_s[3]] = s_list[3, idx_s[2]: idx_s[3]] + correctionshear2
    s56_list[4] = np.absolute(s56_list[1] - correctionshear1 + correctionshear2)

    return s_list


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

    q2 = 1 / 2*(Am2 + Am1*(C2*D1+B2*C1)/(C1*B1+C2*B2))
    q1 = (q2*(C2*D1+B2*C1))/(C1*B1+C2*B2)

    theta_1 = L * C1 * (q1*B1 - q2*B2)
    theta_2 = L * C2 * (q2*D1 - q1*B2)

    theta = L*C1
    print(shearcentre(a))
    return theta_1, theta_2


print(num_twist())