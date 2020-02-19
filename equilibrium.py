import numpy as np

from geometry import *
from loading import *


def distr_defl_func(q1, q2, L, I, E=71000000000):
    # deflection return positive in positive direction of q
    if q1 != q2:
        return L*L*L*L/(E*I)*(q2 /8 + (q1-q2)*3/72)
    else:
        return L*L*L*L*q2/(E*I*8)


def distr_angle_func(q1, q2, L, I, E=71000000000):
    if q1 != q2:
        return L*L*L/(E*I)*(q2/6 + (q1-q2)/18)
    else:
        return L*L*L*q2/(6*E*I)


if __name__ == "__main__":
    a = Dataset()
    Nx = 41
    Nz = 81


    coordlist = make_sections(41, 81, a)
    hinge1_idx = np.argmin(np.abs(coordlist - a.hinge1))
    hinge1_val = coordlist[hinge1_idx]
    hinge2_idx = np.argmin(np.abs(coordlist - a.hinge2))
    hinge2_val = coordlist[hinge2_idx]
    hinge3_idx = np.argmin(np.abs(coordlist - a.hinge3))
    hinge3_val = coordlist[hinge3_idx]

    aerogrid = aero_points(Nx, Nz, a)
    aeroloc, aeroforce = get_aero_resultants("aerodata.csv", aerogrid)
    q = partial(getq, aerogrid, aeroforce, a)
    dx_list = []
    for i in range(coordlist.size - 1):
        dx_list.append(coordlist[i + 1] - coordlist[i])

    u = np.zeros((coordlist.size, 2))
    thetas = np.zeros(coordlist.size)
    for i in range(hinge2_idx, -1, -1):
        u[i, 0] += u[i + 1, 0] + distr_defl_func(-q(coordlist[i]), -q(coordlist[i+1]), dx_list[i], a.Izz())
        thetas[i] = distr_angle_func(-q(coordlist[i]), -q(coordlist[i+1]), dx_list[i], a.Izz())
        for j in range(hinge2_idx-i):
    # TODO: continue here @Max


