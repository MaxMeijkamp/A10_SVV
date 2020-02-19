import numpy as np
import geometry

def von_mises(sigma_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_xz):
    #Takes arg. sigma_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_xz.
    sigma_vm = np.sqrt( ( (sigma_xx - sigma_yy)**2 + (sigma_yy - sigma_zz)**2 + (sigma_zz - sigma_xx)**2 ) / 2 + 3*((tau_xy)**2 + (tau_yz)**2 + (tau_xz)**2) )
    return sigma_vm

def bending(y, z, Iyy, Izz, M_y, M_z, centroid):
    #returns sigma_x for location y, z on the aileron. Inputs: the point x, y;
    #Cont'd: the MOI Izz, Iyy; the relevant moments M_z, M_y (tension in + plane positive), centroid: tuple with xyz centroid of the aileron
    return M_y *  (z - centroid[2]) / Iyy + M_z * (y - centroid[1]) / Izz

def shear():
    return
