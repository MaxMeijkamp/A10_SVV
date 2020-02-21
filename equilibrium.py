import numpy as np

from InputClasses import *
from loading import *
import scipy.linalg as la


def make_u(loads):
    dof_size = 2 if loads.bending == False else 4
    u = np.zeros((loads.gridsize, dof_size))
    thetas = np.zeros((loads.gridsize, int(dof_size/2)))

    for i in range(loads.hinge2_idx - 1, -1, -1):
        u[i, 0] += u[i + 1, 0] + loads.forward_defl_distr(i)
        thetas[i,0] = loads.forward_angle_distr(i)
        if loads.bending:
            u[i, 2] += u[i + 1, 0] + loads.forward_defl_point(i)
            thetas[i,1] = loads.forward_angle_point(i)
        for j in range(1, loads.hinge2_idx - i):
            u[i, 0] += loads.dx_list[i] * thetas[loads.hinge2_idx - j, 0]
            if loads.bending:
                u[i, 2] += loads.dx_list[i] * thetas[loads.hinge2_idx - j, 1]
        u[i, 1] -= loads.dx_list[i]

    for i in range(loads.hinge2_idx + 1, loads.gridsize-1, 1):
        u[i, 0] += u[i - 1, 0] + loads.backward_defl_distr(i)
        thetas[i, 0] = loads.backward_angle_distr(i)
        if loads.bending:
            u[i, 3] += u[i - 1, 0] + loads.backward_defl_point(i)
            thetas[i,1] = loads.backward_angle_point(i)
        for j in range(1, i - loads.hinge2_idx):
            u[i, 0] += loads.dx_list[i] * thetas[loads.hinge2_idx + j, 0]
            if loads.bending:
                u[i, 3] += loads.dx_list[i] * thetas[loads.hinge2_idx + j, 1]
        u[i, 1] += loads.dx_list[i - 1]
    return u



if __name__ == "__main__":
    a = Aileron()
    loads = AppliedLoads(int(a.span*1000), aileron=a) # Section at every millimeter

    u = make_u(loads)
    print(u)
    u1 = u[loads.hinge1_idx]
    u3 = u[loads.hinge3_idx]

    matrix = np.zeros((4,4))
    matrix[0, :-1] = u1[1:]
    matrix[1, :-1] = u3[1:]
    matrix[2, 1:] = np.ones(3)
    matrix[3, 1:-1] = np.array([loads.hinge1_val-loads.hinge2_val, loads.hinge3_val-loads.hinge2_val])
    b = np.array([-u1[0], -u3[0], 0, 0])

    x = la.solve(matrix, b)
    u_new = u[:,0] + u[:,1:].dot(x[:-1])
    print(u_new)
    plt.plot(loads._grid, u_new)
    plt.show()


