import numpy as np
import numericaltools
import InputClasses


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
