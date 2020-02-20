import numpy as np

from InputClasses import *
from loading import *


def make_u(loads):
    dof_size = 2 if loads.bending == False else 4
    u = np.zeros((loads.gridsize, dof_size))
    thetas = np.zeros(loads.gridsize)

    for i in range(loads.hinge_idx(2) - 1, -1, -1):
        u[i, 0] += u[i + 1, 0] + loads.forward_defl_distr(i)
        thetas[i] = loads.forward_angle_distr(i)
        for j in range(1, loads.hinge_idx(2) - i):
            u[i, 0] -= loads.dx_list[i] * thetas[loads.hinge_idx(2) - j]
        u[i, 1] -= loads.dx_list[i]
        if dof_size == 4:
            u[i, 2] += u[i + 1, 0] + loads.forward_defl_point(i)
    for i in range(loads.hinge_idx(2) + 1, loads.gridsize, 1):
        u[i, 0] += u[i - 1, 0] + loads.backward_defl_distr(i)
        thetas[i] = loads.backward_angle_distr(i)
        for j in range(1, i - loads.hinge_idx(2) - i):
            u[i, 0] += loads.dx_list[i] * thetas[loads.hinge_idx(2) + j]
        u[i, 1] += loads.dx_list[i - 1]

    return u



if __name__ == "__main__":
    a = Aileron()
    loads = AppliedLoads(int(a.span*1000), aileron=a) # Section at every millimeter

    u = make_u(a, loads)

    print(u[:15,:])


