import numpy as np
import numericaltools
import InputClasses
from numericaltools import *
import matplotlib.pyplot as plt

def von_mises(sigma_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_xz):
    # Takes arg. sigma_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_xz.
    sigma_xy_dif = sigma_xx - sigma_yy
    sigma_yz_dif = sigma_yy - sigma_zz
    sigma_xz_dif = sigma_zz - sigma_xx
    sigma_dif_tot = (sigma_xy_dif*sigma_xy_dif + sigma_yz_dif*sigma_yz_dif + sigma_xz_dif*sigma_xz_dif) / 2
    tau_tot = 3 * ((tau_xy * tau_xy) + (tau_yz * tau_yz) + (tau_xz * tau_xz))
    sigma_vm = np.sqrt(sigma_dif_tot + tau_tot)
    return sigma_vm


def bending(y, z, Iyy, Izz, M_y, M_z, centroid):
    # returns sigma_x for location y, z on the aileron. Inputs: the point x, y; the MOI Izz, Iyy;
    # the relevant moments M_z, M_y (tension in + plane positive), centroid: tuple with xyz centroid of the aileron
    return M_y * (z - centroid[2]) / Iyy + M_z * (y - centroid[1]) / Izz

# For bending, not only sigma_x is needed, but also the displacements and angular displacements
# at all sections caused by the bending.

def shearflow(dataset, shearforce = 1.):
    # Calculates shear flows around circumference and in spar due to shear force applied in SC
    aileron = dataset
    r = aileron.radius
    a = aileron.a
    Izz = aileron.Izz()
    skint = aileron.skint
    spart = aileron.spart

    dist_s1 = dist_s4 = .5 * np.pi * r
    dist_s2 = dist_s3 = a
    dist_s5 = dist_s6 = r

    def func_s(s):
        if s < dist_s1:
            return (r * np.sin(s / r))
        elif dist_s1 <= s < dist_s1 + dist_s2:
            return (r - (s - dist_s1) * r / a)
        elif dist_s2 + dist_s1 <= s < dist_s3 + dist_s2 + dist_s1:
            return (0. - (s - dist_s1 - dist_s2) * r / a)
        elif dist_s3 + dist_s2 + dist_s1 <= s <= dist_s4 + dist_s3 + dist_s2 + dist_s1:
            return - r * np.cos((s - dist_s1 - dist_s2 - dist_s3) / r)

    def func_s56(s):  # from s = 0 to s = height / 2
        return s

    integrateV = np.vectorize(integrate)
    func_s = np.vectorize(func_s)
    func_s56 = np.vectorize(func_s56)
    num_steps = 2000

    s_list = np.zeros((7, num_steps))  # array with s value, integral shear, booms
    s_list[0] = np.linspace(0, aileron._circumference, num_steps)
    idx_s = [
        0, int((np.abs(s_list[0] - dist_s1)).argmin()),
        int(np.abs(s_list[0] - dist_s1 - dist_s2).argmin()),
        int(np.abs(s_list[0] - dist_s1 - dist_s2 - dist_s3).argmin()),
        int(np.abs(s_list[0] - dist_s1 - dist_s2 - dist_s3 - dist_s4).argmin())
    ]

    ss = np.arange(-r, r, aileron._circumference / num_steps)
    s56_list = np.zeros((7, ss.size))
    s56_list[0] = ss

    # Skin base shear flows
    s56_list[1, 0:int(ss.size / 2)] = 1. * shearforce / Izz * (
        integrateV(func_s56, 0, s56_list[0, 0:int(ss.size / 2)], 2000)) * spart
    s56_list[1, int(ss.size / 2):int(ss.size)] = -1. * shearforce / Izz * (
        integrateV(func_s56, 0, s56_list[0, int(ss.size / 2):int(ss.size)], 2000)) * spart
    s_list[1] = -1. / Izz * shearforce * (integrateV(func_s, 0, s_list[0], 5000)) * skint

    # Base shear flow of booms:
    boom_stepsize = aileron._circumference / aileron.stiffn
    for i in range(aileron.stiffn):
        boomdist = boom_stepsize * i
        idx = (np.abs(s_list[0] - boomdist)).argmin()
        s_list[2, idx: idx_s[4]] += (-1. * shearforce / Izz * aileron.stiffener_area * func_s(s_list[0, idx]))

    # Booms and skin contribution combined in row 3 of s_list:
    s_list[3] = s_list[2] + s_list[1]

    # Find correction shear flow:
    stepwidth = aileron._circumference / num_steps  # s_list[0,1] - s_list[0,0]

    shearbasesum1 = (stepwidth * ((np.sum(s_list[3, 0: idx_s[1]]) +
                                   np.sum(s_list[3, idx_s[3]: idx_s[4]])) / skint -
                                  np.sum(s56_list[1]) / spart))
    shearbasesum2 = (stepwidth * ((np.sum(s_list[3, idx_s[1]: idx_s[3]])) / skint +
                                  np.sum(s56_list[1]) / spart))

    twistmatrix = np.array([[2 * dist_s1 / skint + 2 * dist_s5 / spart, -2 * dist_s5 / spart],
                            [-2 * dist_s5 / spart, 2 * dist_s2 / skint + 2 * dist_s5 / spart]])

    shearbasevector = np.array([[-shearbasesum1],
                                [-shearbasesum2]])
    correctionshears = np.linalg.solve(twistmatrix, shearbasevector)
    correctionshear1 = correctionshears[0]
    correctionshear2 = correctionshears[1]

    # check should be 0 vector
    check = twistmatrix @ correctionshears - shearbasevector

    # Add correction shears and base shear to row 4
    s_list[4, idx_s[0]: idx_s[1]] = s_list[3, idx_s[0]: idx_s[1]] + correctionshear1
    s_list[4, idx_s[3]: idx_s[4]] = s_list[3, idx_s[3]: idx_s[4]] + correctionshear1
    s_list[4, idx_s[1]: idx_s[3]] = s_list[3, idx_s[1]: idx_s[3]] + correctionshear2
    s56_list[4] = s56_list[1] - correctionshear1 + correctionshear2

    return s_list[[0,4,3]], s56_list[[0,4,1]]

if __name__ == '__main__':
    flow = shearflow(InputClasses.Aileron())[0]
    plt.plot(flow[0], flow[1])
    plt.show()