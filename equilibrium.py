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
        # Previous known deflection plus aero distributed load
        u[i, 0] += u[i + 1, 0] + loads.forward_defl_distr(i)
        # Contribution of internal shear and moment due to aero load
        u[i, 0] += -loads.pointload_defl(shearforces[i,0], loads.dx_list[i], loads.a.Izz()) + loads.moment_defl(moments[i,0], loads.dx_list[i], loads.a.Izz())
        # Angle due to aero load
        thetas[i,0] = loads.forward_angle_distr(i)
        # Angle due to internal shear and moment due to aero load
        thetas[i,0] += loads.moment_angle(moments[i,0], loads.dx_list[i], loads.a.Izz()) - loads.pointload_angle(shearforces[i,0], loads.dx_list[i], loads.a.Izz())

        # Previous deflection due to angle at D and new contribution
        u[i, 1] += u[i+1,1] - loads.dx_list[i]

        if loads.bending:
            # Previous deflection plus contribution of force 1
            u[i, 2] += u[i + 1, 2] + loads.forward_defl_point(i)
            # Adding deflection due to internal shear and moment due to force 1
            u[i, 2] += -loads.pointload_defl(shearforces[i,1], loads.dx_list[i], loads.a.Izz()) + loads.moment_defl(moments[i,1], loads.dx_list[i], loads.a.Izz())
            # Angle due to force 1
            thetas[i,1] = loads.forward_angle_point(i)
            # Angle due to internal shear and moment due to force 1
            thetas[i,1] += -loads.pointload_angle(shearforces[i,1], loads.dx_list[i], loads.a.Izz()) + loads.moment_angle(moments[i,1], loads.dx_list[i], loads.a.Izz())

        for j in range(1, loads.hinge2_idx - i):
            # Add contribution of the angle depending on the aero load
            u[i, 0] += loads.dx_list[i] * thetas[loads.hinge2_idx - j, 0]
            if loads.bending:
                # Add contribution of the angle depending on force 1
                u[i, 2] += loads.dx_list[i] * thetas[loads.hinge2_idx - j, 1]

    for i in range(loads.hinge2_idx + 1, loads.gridsize, 1):
        # Previous known deflection plus aero distributed load
        u[i, 0] += u[i - 1, 0] + loads.backward_defl_distr(i)
        # Contribution of internal shear and moment due to aero load
        u[i, 0] += -loads.pointload_defl(shearforces[i-1,0], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_defl(moments[i-1,0], loads.dx_list[i-1], loads.a.Izz())
        # Angle due to aero load
        thetas[i, 0] = loads.backward_angle_distr(i)
        # Angle due to internal shear and moment due to aero load
        thetas[i, 0] += loads.moment_angle(moments[i-1,0], loads.dx_list[i-1], loads.a.Izz()) - loads.pointload_angle(shearforces[i-1,0], loads.dx_list[i-1], loads.a.Izz())

        # Previous deflection due to angle at D and new contribution
        u[i, 1] += u[i-1,1] + loads.dx_list[i - 1]


        u[i, 3] += u[i-1,3] - loads.pointload_defl(shearforces[i-1,2], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_defl(moments[i-1,2], loads.dx_list[i-1], loads.a.Izz())
        thetas[i,2] = -loads.pointload_angle(shearforces[i-1,2], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_angle(moments[i-1,2], loads.dx_list[i-1], loads.a.Izz())

        if loads.bending:
            u[i, 4] += u[i - 1, 4] + loads.backward_defl_point(i)
            u[i, 4] += -loads.pointload_defl(shearforces[i-1,3], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_defl(moments[i-1,3], loads.dx_list[i-1], loads.a.Izz())
            thetas[i,3] = loads.backward_angle_point(i)
            thetas[i,3] += -loads.pointload_angle(shearforces[i-1,3], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_angle(moments[i-1,3], loads.dx_list[i-1], loads.a.Izz())

            u[i, 2] += u[i - 1, 2]
            u[i, 2] += -loads.pointload_defl(shearforces[i-1,1], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_defl(moments[i-1,1], loads.dx_list[i-1], loads.a.Izz())
            thetas[i,1] = -loads.pointload_angle(shearforces[i-1,1], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_angle(moments[i-1,1], loads.dx_list[i-1], loads.a.Izz())

        for j in range(1, i - loads.hinge2_idx):
            u[i, 0] += loads.dx_list[i-1] * thetas[loads.hinge2_idx + j, 0]
            u[i, 3] += loads.dx_list[i-1] * thetas[loads.hinge2_idx + j, 2]
            if loads.bending:
                u[i, 2] += loads.dx_list[i-1] * thetas[loads.hinge2_idx + j, 1]
                u[i, 4] += loads.dx_list[i-1] * thetas[loads.hinge2_idx + j, 3]

    return u, shearforces, moments


def make_u_z(loads):
    u = np.zeros((loads.gridsize, 6))
    thetas = np.zeros((loads.gridsize, 5))
    shearforces = loads.int_shear_z(loads.grid).T
    moments = loads.int_moment_z(loads.grid).T
    for i in range(loads.hinge2_idx - 1, -1, -1):
        # Previous known deflection plus internal shear due to P
        u[i, 0] += u[i + 1, 0] + loads.pointload_defl(shearforces[i,0], loads.dx_list[i], loads.a.Iyy())
        # Contribution of internal moment due to P, negative due to positive definition of internal moment
        u[i, 0] -= loads.moment_defl(moments[i,0], loads.dx_list[i], loads.a.Iyy())
        # Angle due to internal shear and moment due to P
        thetas[i,0] = loads.pointload_angle(shearforces[i,0], loads.dx_list[i], loads.a.Iyy()) - loads.moment_angle(moments[i,0], loads.dx_list[i], loads.a.Iyy())

        # Previous deflection due to angle at D and new contribution
        u[i, 1] += u[i+1,1] - loads.dx_list[i]

        if loads.bending:
            # Previous deflection plus contribution of force 1
            u[i, 2] += u[i + 1, 2] + loads.defl_point_z(i, 0, 0)
            # Adding deflection due to internal shear and moment due to force 1
            u[i, 2] += loads.pointload_defl(shearforces[i,1], loads.dx_list[i], loads.a.Iyy()) - loads.moment_defl(moments[i,1], loads.dx_list[i], loads.a.Iyy())
            # Angle due to force 1
            thetas[i,1] = loads.angle_point_z(i, 0, 0)
            # Angle due to internal shear and moment due to force 1
            thetas[i,1] += loads.pointload_angle(shearforces[i,1], loads.dx_list[i], loads.a.Iyy()) - loads.moment_angle(moments[i,1], loads.dx_list[i], loads.a.Iyy())
        if loads.jammed:
            # Previous deflection plus contribution of force a
            u[i, 3] += u[i + 1, 3] + loads.defl_point_z(i, 1, 0)
            # Adding deflection due to internal shear and moment due to force a
            u[i, 3] += loads.pointload_defl(shearforces[i, 2], loads.dx_list[i], loads.a.Iyy()) - loads.moment_defl(
                moments[i, 2], loads.dx_list[i], loads.a.Iyy())
            # Angle due to force a
            thetas[i, 2] = loads.angle_point_z(i, 1, 0)
            # Angle due to internal shear and moment due to force a
            thetas[i, 2] += loads.pointload_angle(shearforces[i, 2], loads.dx_list[i], loads.a.Iyy()) - loads.moment_angle(moments[i, 2], loads.dx_list[i], loads.a.Iyy())

        for j in range(1, loads.hinge2_idx - i):
            # Add contribution of the angle depending on P
            u[i, 0] += loads.dx_list[i] * thetas[loads.hinge2_idx - j, 0]
            if loads.bending:
                # Add contribution of the angle depending on force 1
                u[i, 2] += loads.dx_list[i] * thetas[loads.hinge2_idx - j, 1]
            if loads.jammed:
                # Add contribution of the angle depending on force 1
                u[i, 3] += loads.dx_list[i] * thetas[loads.hinge2_idx - j, 2]

    for i in range(loads.hinge2_idx + 1, loads.gridsize, 1):
        # Previous known deflection plus internal shear due to P
        u[i, 0] += u[i - 1, 0] + loads.pointload_defl(shearforces[i-1,0], loads.dx_list[i-1], loads.a.Iyy())
        # Contribution of internal shear and moment due to aero load
        u[i, 0] -= loads.moment_defl(moments[i-1,0], loads.dx_list[i-1], loads.a.Iyy())
        # Angle due to internal shear and moment due to P
        thetas[i, 0] = loads.pointload_angle(shearforces[i-1,0], loads.dx_list[i-1], loads.a.Iyy()) - loads.moment_angle(moments[i-1,0], loads.dx_list[i-1], loads.a.Iyy())

        # Previous deflection due to angle at D and new contribution
        u[i, 1] += u[i-1,1] + loads.dx_list[i - 1]


        u[i, 4] += u[i-1,4] + loads.pointload_defl(shearforces[i-1,3], loads.dx_list[i-1], loads.a.Iyy()) - loads.moment_defl(moments[i-1,3], loads.dx_list[i-1], loads.a.Iyy())
        thetas[i,3] = loads.pointload_angle(shearforces[i-1,3], loads.dx_list[i-1], loads.a.Iyy()) - loads.moment_angle(moments[i-1,3], loads.dx_list[i-1], loads.a.Iyy())

        if loads.bending:
            u[i, 5] += u[i - 1, 5] + loads.defl_point_z(i, 4, 1)
            u[i, 5] += loads.pointload_defl(shearforces[i-1,4], loads.dx_list[i-1], loads.a.Iyy()) - loads.moment_defl(moments[i-1,4], loads.dx_list[i-1], loads.a.Iyy())
            thetas[i,4] = loads.angle_point_z(i, 4, 1)
            thetas[i,4] += loads.pointload_angle(shearforces[i-1,4], loads.dx_list[i-1], loads.a.Iyy()) - loads.moment_angle(moments[i-1,4], loads.dx_list[i-1], loads.a.Iyy())

            u[i, 2] += u[i - 1, 2]
            u[i, 2] += loads.pointload_defl(shearforces[i-1,1], loads.dx_list[i-1], loads.a.Iyy()) - loads.moment_defl(moments[i-1,1], loads.dx_list[i-1], loads.a.Iyy())
            thetas[i,1] = loads.pointload_angle(shearforces[i-1,1], loads.dx_list[i-1], loads.a.Iyy()) - loads.moment_angle(moments[i-1,1], loads.dx_list[i-1], loads.a.Iyy())
        if loads.jammed:
            # Previous deflection plus contribution of force a
            u[i, 3] += u[i - 1, 3] + loads.defl_point_z(i, 1, 1)
            # Adding deflection due to internal shear and moment due to force a
            u[i, 3] += loads.pointload_defl(shearforces[i-1, 2], loads.dx_list[i-1], loads.a.Iyy()) - loads.moment_defl(moments[i-1, 2], loads.dx_list[i-1], loads.a.Iyy())
            # Angle due to force a
            thetas[i, 2] = loads.angle_point_z(i, 1, 1)
            # Angle due to internal shear and moment due to force a
            thetas[i, 2] += loads.pointload_angle(shearforces[i-1, 2], loads.dx_list[i-1], loads.a.Iyy()) - loads.moment_angle(moments[i-1, 2], loads.dx_list[i-1], loads.a.Iyy())


        for j in range(1, i - loads.hinge2_idx):
            u[i, 0] += loads.dx_list[i-1] * thetas[loads.hinge2_idx + j, 0]
            u[i, 4] += loads.dx_list[i-1] * thetas[loads.hinge2_idx + j, 3]
            if loads.bending:
                u[i, 2] += loads.dx_list[i-1] * thetas[loads.hinge2_idx + j, 1]
                u[i, 5] += loads.dx_list[i-1] * thetas[loads.hinge2_idx + j, 4]
            if loads.jammed:
                # Add contribution of the angle depending on force 1
                u[i, 3] += loads.dx_list[i-1] * thetas[loads.hinge2_idx - j, 2]

    return u, shearforces, moments


def calc_deflection_y(loads, u=np.nan):
    if np.any(np.isnan(u)):
        u = make_u_y(loads)
    u1 = u[loads.hinge1_idx]
    u3 = u[loads.hinge3_idx]

    matrix = np.zeros((4,4))
    matrix[0, :] = u1[1:]
    matrix[1, :] = u3[1:]
    matrix[2, 1:] = np.ones(3)
    matrix[3, 1:] = np.array([loads.hinge1_val-loads.hinge2_val, 0, loads.hinge3_val-loads.hinge2_val])
    b = np.array([loads.hinge1d_y - u1[0], loads.hinge3d_y - u3[0], loads.totforce,
                  loads.totforce * (loads.res_xloc - loads.hinge2_val)])

    x = la.solve(matrix, b) # Theta_d, F_y1, F_y2, F_y3
    u_solved = u[:,0] + u[:,1:].dot(x)
    return u_solved, x


def calc_deflection_z(loads, u=np.nan):
    if np.any(np.isnan(u)):
        u = make_u_y(loads)
    u1 = u[loads.hinge1_idx]
    ua = u[loads.act1_idx]
    u3 = u[loads.hinge3_idx]

    matrix = np.zeros((5,5))
    matrix[0, :] = u1[1:]
    matrix[1, :] = ua[1:]
    matrix[2, :] = u3[1:]
    matrix[3, 1:] = np.ones(4)
    matrix[4, 1:] = np.array([loads.hinge1_val-loads.hinge2_val, loads.act1_val-loads.hinge2_val, 0, -loads.hinge3_val-loads.hinge2_val])
    b = np.array([- u1[0], - ua[0], - u3[0], loads.P, loads.P * (loads.act2_val-loads.hinge2_val)])

    x = la.solve(matrix, b) # Theta_d, F_y1, F_y_I, F_y2, F_y3
    u_solved = u[:,0] + u[:,1:].dot(x)
    return u_solved, x


def alt_u_y(loads):
    u = np.zeros((loads.gridsize, 6))
    thetas = np.zeros((loads.gridsize, 4))
    shearforces = loads.int_shear_y(loads.grid).T
    moments = loads.int_moment_y(loads.grid).T

    u[0, :-1] = np.zeros(u[0,:-1].size)
    u[:, -1] = np.ones(u[:, -1].size)
    for i in range(1, loads.gridsize, 1):
        # Previous known deflection plus aero distributed load
        u[i, 0] += u[i - 1, 0] + loads.backward_defl_distr(i)
        # Contribution of internal shear and moment due to aero load
        u[i, 0] += loads.pointload_defl(shearforces[i-1,0], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_defl(moments[i-1,0], loads.dx_list[i-1], loads.a.Izz())
        # Angle due to aero load
        thetas[i, 0] = loads.backward_angle_distr(i)
        # Angle due to internal shear and moment due to aero load
        thetas[i, 0] += loads.moment_angle(moments[i-1,0], loads.dx_list[i-1], loads.a.Izz()) + loads.pointload_angle(shearforces[i-1,0], loads.dx_list[i-1], loads.a.Izz())

        # Previous deflection due to angle at D and new contribution
        u[i, 1] += u[i-1,1] + loads.dx_list[i - 1]


        u[i, 3] += u[i-1,3] + loads.pointload_defl(shearforces[i-1,2], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_defl(moments[i-1,2], loads.dx_list[i-1], loads.a.Izz())
        thetas[i,2] = loads.pointload_angle(shearforces[i-1,2], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_angle(moments[i-1,2], loads.dx_list[i-1], loads.a.Izz())

        if loads.bending:
            u[i, 4] += u[i - 1, 4] + loads.backward_defl_point(i)
            u[i, 4] += loads.pointload_defl(shearforces[i-1,3], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_defl(moments[i-1,3], loads.dx_list[i-1], loads.a.Izz())
            thetas[i,3] = loads.backward_angle_point(i)
            thetas[i,3] += loads.pointload_angle(shearforces[i-1,3], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_angle(moments[i-1,3], loads.dx_list[i-1], loads.a.Izz())

            u[i, 2] += u[i - 1, 2]
            u[i, 2] += loads.pointload_defl(shearforces[i-1,1], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_defl(moments[i-1,1], loads.dx_list[i-1], loads.a.Izz())
            thetas[i,1] = loads.pointload_angle(shearforces[i-1,1], loads.dx_list[i-1], loads.a.Izz()) + loads.moment_angle(moments[i-1,1], loads.dx_list[i-1], loads.a.Izz())

        for j in range(1, i):
            u[i, 0] += loads.dx_list[i-1] * thetas[j, 0]
            u[i, 3] += loads.dx_list[i-1] * thetas[j, 2]
            if loads.bending:
                u[i, 2] += loads.dx_list[i-1] * thetas[j, 1]
                u[i, 4] += loads.dx_list[i-1] * thetas[j, 3]

    return u, shearforces, moments


def alt_y_defl(loads, u=np.nan):
    if np.any(np.isnan(u)):
        u = alt_u_y(loads)
    u1 = u[loads.hinge1_idx]
    u2 = u[loads.hinge2_idx]
    u3 = u[loads.hinge3_idx]

    matrix = np.zeros((5,5))
    matrix[0, :] = u1[1:]
    matrix[1, :] = u2[1:]
    matrix[2, :] = u3[1:]
    matrix[3, 1:-1] = np.ones(3)
    matrix[4, 1:-1] = np.array([loads.hinge1_val-loads.hinge2_val, 0, loads.hinge3_val-loads.hinge2_val])
    b = np.array([loads.hinge1d_y - u1[0], -u2[0], loads.hinge3d_y - u3[0], loads.totforce,
                  loads.totforce * (loads.res_xloc - loads.hinge2_val)])

    x = la.solve(matrix, b) # Theta_d, F_y1, F_y2, F_y3, u_0_0
    u_solved = u[:,0] + u[:,1:].dot(x)
    return u_solved, x



def plot_u(loads, u, axis, equal=False) -> None:
    # axis = 0 is y axis, axis = 1 is z axis
    if axis == 0:
        plt.plot(loads.grid, u)
        plt.scatter(loads.hinge1_val, loads.hinge1d_y)
        plt.scatter(loads.hinge2_val, 0)
        plt.scatter(loads.hinge3_val, loads.hinge3d_y)
    elif axis == 1:
        plt.plot(loads.grid, u)
        plt.scatter(loads.hinge1_val, 0)
        plt.scatter(loads.act1_val, 0)
        plt.scatter(loads.hinge2_val, 0)
        plt.scatter(loads.hinge3_val, 0)
    if equal:
        plt.axis('equal')
    plt.show()

def save_data(x, appendix = "", prefix = "", shear_blank=False, shear_solved=False, moment_blank=False, moment_solved=False, u_blank=False, u=False) -> None:
    if appendix != "":
        appendix = "_"+appendix
    np.savetxt(prefix+"x"+appendix+".csv", x, delimiter=",")
    if np.any(shear_blank) != False:
        np.savetxt(prefix+"shearblank"+appendix+".csv", shear_blank, delimiter=",")
    if np.any(shear_solved) != False:
        np.savetxt(prefix+"shearsolved"+appendix+".csv", shear_solved, delimiter=",")
    if np.any(moment_blank) != False:
        np.savetxt(prefix+"momentblank"+appendix+".csv", moment_blank, delimiter=",")
    if np.any(moment_solved) != False:
        np.savetxt(prefix+"momentsolved"+appendix+".csv", moment_solved, delimiter=",")
    if np.any(u_blank) != False:
        np.savetxt(prefix+"ublank"+appendix+".csv", u_blank, delimiter=",")
    if np.any(u) != False:
        np.savetxt(prefix+"usolved"+appendix+".csv", u, delimiter=",")

if __name__ == "__main__":
    a = Aileron()
    loads = AppliedLoads(int(a.span*1000), aileron=a) # Section at every millimeter

    u, shear, moment = make_u_z(loads)
    # np.savetxt("shear2.csv", shear, delimiter=",")
    # np.savetxt("moment2.csv", moment, delimiter=",")


    # print(u[[x+loads.hinge3_idx for x in range(-3,4)]])
    u_new, x = calc_deflection_z(loads, u)

    shearsolved = shear[:,0] + shear[:, 1:].dot(x[1:])
    momentsolved = shear[:,0] + shear[:, 1:].dot(x[1:])
    # np.savetxt("shearsolved.csv", shearsolved, delimiter=",")

    print(x)
    u_save = np.copy(u)
    u_save[:,1] *= x[0]
    u_save[:,2] *= x[1]
    u_save[:,3] *= x[2]
    u_save[:,4] *= x[3]
    # u_save[:,5] *= x[4]
    # np.savetxt("u_new_data3.csv", u_save, delimiter=",")
    # u_y_new, x = calc_deflection_y(loads)
    #
    # print(u_y_new[loads.hinge2_idx])
    #
    plot_u(loads, u_new, 1)


