import numpy as np

from InputClasses import *
from loading import *
import scipy.linalg as la


def make_u_y(loads):
    u = np.zeros((loads.gridsize, 5))
    thetas = np.zeros((loads.gridsize, 4))
    shearforces = loads.int_shear_y(loads.grid).T
    moments = loads.int_moment_y(loads.grid).T
    for i in range(loads.hinge2_idx - 1, -1, -1):
        u[i, 0] += u[i + 1, 0] + loads.forward_defl_distr(i)
        u[i, 0] += loads.pointload_defl(shearforces[i,0], loads.dx_list[i], loads.a.Izz()) + loads.moment_defl(moments[i,0], loads.dx_list[i], loads.a.Izz())
        thetas[i,0] = loads.forward_angle_distr(i)
        thetas[i,0] += loads.moment_angle(shearforces[i,0], loads.dx_list[i], loads.a.Izz()) + loads.pointload_angle(shearforces[i,0], loads.dx_list[i], loads.a.Izz())

        u[i, 1] += u[i+1,1] - loads.dx_list[i]

        if loads.bending:
            u[i, 2] += u[i + 1, 2] + loads.forward_defl_point(i)
            u[i, 2] += loads.pointload_defl(shearforces[i,1], loads.dx_list[i], loads.a.Izz()) + loads.moment_defl(moments[i,1], loads.dx_list[i], loads.a.Izz())
            thetas[i,1] = loads.forward_angle_point(i)
            thetas[i,1] += loads.pointload_angle(shearforces[i,1], loads.dx_list[i], loads.a.Izz()) + loads.moment_angle(moments[i,1], loads.dx_list[i], loads.a.Izz())

        for j in range(1, loads.hinge2_idx - i):
            u[i, 0] += loads.dx_list[i] * thetas[loads.hinge2_idx - j, 0]
            if loads.bending:
                u[i, 2] += loads.dx_list[i] * thetas[loads.hinge2_idx - j, 1]

    for i in range(loads.hinge2_idx + 1, loads.gridsize, 1):
        u[i, 0] += u[i - 1, 0] + loads.backward_defl_distr(i)
        u[i, 0] += loads.pointload_defl(shearforces[i,0], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_defl(moments[i,0], loads.dx_list[i-1], loads.a.Izz())
        thetas[i, 0] = loads.backward_angle_distr(i)
        thetas[i, 0] += loads.moment_angle(shearforces[i,0], loads.dx_list[i-1], loads.a.Izz()) + loads.pointload_angle(shearforces[i,0], loads.dx_list[i-1], loads.a.Izz())

        u[i, 1] += u[i-1,1] + loads.dx_list[i - 1]

        u[i, 3] += u[i-1,3] + loads.pointload_defl(shearforces[i,2], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_defl(moments[i,2], loads.dx_list[i-1], loads.a.Izz())
        thetas[i,2] = loads.pointload_angle(shearforces[i,2], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_angle(moments[i,2], loads.dx_list[i-1], loads.a.Izz())

        if loads.bending:
            u[i, 4] += u[i - 1, 4] + loads.backward_defl_point(i)
            u[i, 4] += loads.pointload_defl(shearforces[i,3], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_defl(moments[i,3], loads.dx_list[i-1], loads.a.Izz())
            thetas[i,3] = loads.backward_angle_point(i)
            thetas[i,3] += loads.pointload_angle(shearforces[i,3], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_angle(moments[i,3], loads.dx_list[i-1], loads.a.Izz())

            u[i, 2] += u[i - 1, 2]
            u[i, 2] += loads.pointload_defl(shearforces[i,1], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_defl(moments[i,1], loads.dx_list[i-1], loads.a.Izz())
            thetas[i,1] = loads.pointload_angle(shearforces[i,1], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_angle(moments[i,1], loads.dx_list[i-1], loads.a.Izz())

        for j in range(1, i - loads.hinge2_idx):
            u[i, 0] += loads.dx_list[i-1] * thetas[loads.hinge2_idx + j, 0]
            u[i, 3] += loads.dx_list[i-1] * thetas[loads.hinge2_idx + j, 2]
            if loads.bending:
                u[i, 2] += loads.dx_list[i-1] * thetas[loads.hinge2_idx - j, 1]
                u[i, 4] += loads.dx_list[i-1] * thetas[loads.hinge2_idx + j, 3]

    return u, shearforces, moments

def make_u_z(loads):
    dof_size = 6 if loads.jammed else 5
    u = np.zeros((loads.gridsize, dof_size))
    thetas = np.zeros((loads.gridsize, 4 if loads.jammed else 3))

    for i in range(loads.hinge2_idx - 1, -1, -1):
        u[i, 0] += u[i + 1, 0] + loads.forward_defl_distr(i)
        thetas[i,0] = loads.forward_angle_distr(i)
        if loads.jammed:
            u[i, 2] += u[i + 1, 0] + loads.forward_defl_point(i)
            thetas[i,1] = loads.forward_angle_point(i)
        for j in range(1, loads.hinge2_idx - i):
            u[i, 0] += loads.dx_list[i] * thetas[loads.hinge2_idx - j, 0]
            if loads.jammed:
                u[i, 2] += loads.dx_list[i] * thetas[loads.hinge2_idx - j, 1]
        u[i, 1] += u[i+1,1] - loads.dx_list[i]

    for i in range(loads.hinge2_idx + 1, loads.gridsize, 1):
        u[i, 0] += u[i - 1, 0] + loads.backward_defl_distr(i)
        thetas[i, 0] = loads.backward_angle_distr(i)
        if loads.jammed:
            u[i, 3] += u[i - 1, 0] + loads.backward_defl_point(i)
            thetas[i,1] = loads.backward_angle_point(i)
        for j in range(1, i - loads.hinge2_idx):
            u[i, 0] += loads.dx_list[i-1] * thetas[loads.hinge2_idx + j, 0]
            if loads.jammed:
                u[i, 3] += loads.dx_list[i-1] * thetas[loads.hinge2_idx + j, 1]
        u[i, 1] += u[i-1,1] + loads.dx_list[i - 1]
    return u

def calc_deflection_y(loads, u=np.nan):
    if np.any(np.isnan(u)):
        u = make_u_y(loads)
    u1 = u[loads.hinge1_idx]
    u3 = u[loads.hinge3_idx]

    matrix = np.zeros((4,4))
    matrix[0, :] = u1[1:]
    matrix[1, :] = u3[1:]
    matrix[2, 1:] = np.ones(3)
    matrix[3, 1:-1] = np.array([loads.hinge1_val-loads.hinge2_val, loads.hinge3_val-loads.hinge2_val])
    b = np.array([loads.hinge1d_y - u1[0], loads.hinge3d_y - u3[0], loads.totforce,
                  loads.totforce * (loads.res_xloc - loads.hinge2_val)])

    x = la.solve(matrix, b) # Theta_d, F_y1, F_y2, F_y3
    u_solved = u[:,0] + u[:,1:].dot(x)
    return u_solved, x

def calc_deflection_z(loads):
    u = make_u_z(loads)
    u1 = u[loads.hinge1_idx]
    u3 = u[loads.hinge3_idx]
    uI = u[loads.act1_idx]

    matrix = np.zeros((5,5))
    matrix[0, :] = u1[1:]
    matrix[1, :] = u3[1:]
    matrix[2, 1:] = np.ones(3)
    matrix[3, 1:-1] = np.array([loads.hinge1_val-loads.hinge2_val, loads.hinge3_val-loads.hinge2_val])
    b = np.array([loads.hinge1d_y-u1[0], loads.hinge3d_y-u3[0], loads.totforce, loads.totforce * (loads.res_xloc - loads.hinge2_val)])

    x = la.solve(matrix, b) # Theta_d, F_y1, F_y3, F_y2
    u_solved = u[:,0] + u[:,1:].dot(x[:-1])
    return u_solved, x



if __name__ == "__main__":
    a = Aileron()
    loads = AppliedLoads(int(a.span*1000), aileron=a) # Section at every millimeter

    u, shear, moment = make_u_y(loads)
    # np.savetxt("shear.csv", shear, delimiter=",")
    # np.savetxt("moment.csv", moment, delimiter=",")

    # print(u[[x+loads.hinge3_idx for x in range(-3,4)]])
    u_new, x = calc_deflection_y(loads, u)
    print(x)
    u_save = np.copy(u)
    u_save[:,1] *= x[0]
    u_save[:,2] *= x[1]
    u_save[:,3] *= x[2]
    u_save[:,4] *= x[3]
    # np.savetxt("u_new_data.csv", u_save, delimiter=",")
    # u_y_new, x = calc_deflection_y(loads)
    #
    # print(u_y_new[loads.hinge2_idx])
    #
    plt.plot(loads.grid, u_new)
    plt.scatter(loads.hinge1_val, loads.hinge1d_y)
    plt.scatter(loads.hinge2_val, 0)
    plt.scatter(loads.hinge3_val, loads.hinge3d_y)
    # plt.axis('equal')
    plt.show()


