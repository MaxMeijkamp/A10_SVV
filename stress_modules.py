import numpy as np
import numericaltools
import InputClasses


def von_mises(sigma_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_xz):
    # Takes arg. sigma_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_xz.
    sigma_vm = np.sqrt( ( (sigma_xx - sigma_yy)**2 + (sigma_yy - sigma_zz)**2 + (sigma_zz - sigma_xx)**2 ) / 2 + 3*((tau_xy)**2 + (tau_yz)**2 + (tau_xz)**2) )
    return sigma_vm


def bending(y, z, Iyy, Izz, M_y, M_z, centroid):
    # returns sigma_x for location y, z on the aileron. Inputs: the point x, y; the MOI Izz, Iyy;
    # the relevant moments M_z, M_y (tension in + plane positive), centroid: tuple with xyz centroid of the aileron
    return M_y * (z - centroid[2]) / Iyy + M_z * (y - centroid[1]) / Izz


def find_shear_center(dataset):
    Izz = dataset.Izz()
    SC_y = dataset.centroid[1]
    Vy = 1  # force applied in shear centre
    # cuts at front of aileron and middle in spar
    # y,z axis on middle of spar
    r = dataset._radius
    a = dataset.a
    skint = dataset.skint
    height = dataset.maxy
    stiffenerA = dataset._stiffener_area
    # funcs will return y as function of s
    def func_s1(s): # from s = 0 to s = half a quarter
        return (r * np.sin(s / (2* np.pi * r))) * skint
    def func_s4(s): # from s = 0 to s = half a quarter
        return (r * np.sin(s / (2* np.pi * r)) - r) * skint
    def func_s2(s): # from s = 0 to s = dataset.a
        return (dataset.maxy - s * abs(dataset.maxy) / abs(dataset.maxz)) * skint
    def func_s3(s): # from s = 0 to s = dataset.a
        return (0. - s * abs(dataset.maxy) / abs(dataset.maxz)) * skint
    def func_s5(s): # from s = 0 to s = height / 2
        return (dataset.miny + s) * skint
    def func_s6(s): # from s = 0 to s = height / 2
        return s * skint
    dist_s1, dist_s4 = .5 * np.pi * r, .5 * np.pi * r
    dist_s2, dist_s3 = dataset.a, dataset.a
    dist_s5, dist_s6 = dataset.maxy
    # now per section (1 to 6), find the total q_base distribution:
    # after the integration The boom influences are PUT IN MANUALLY, PUT IN MANUALLY, PUT IN MANUALLY!!!

    dqb1 = -1. / Izz * (numericaltools.integrate(func_s1, 0, dist_s1, 100000)     )
    dqb2 = -1. / Izz * (numericaltools.integrate(func_s2, 0, dist_s2, 100000)    )
    dqb3 = -1. / Izz * (numericaltools.integrate(func_s3, 0, dist_s3, 100000)    )
    dqb4 = -1. / Izz * (numericaltools.integrate(func_s4, 0, dist_s4, 100000)    )
    dqb5 = -1. / Izz * (numericaltools.integrate(func_s5, 0, dist_s5, 100000)    )
    dqb6 = -1. / Izz * (numericaltools.integrate(func_s6, 0, dist_s6, 100000)    )

    return
# For bending, not only sigma_x is needed, but also the displacements and angular displacements
# at all sections caused by the bending.

