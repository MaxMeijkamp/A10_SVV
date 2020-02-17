import numpy as np

simga_xx = 1
sigma_yy = 2
sigma_zz = 1
tau_xy = 2
tau_yz = 1
tau_xz = 2

def von_mises(simga_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_xz)
    #Takes arg. simga_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_xz
    sigma_vm = np.sqrt( ( (sigma_xx - sigma_yy)**2 + (sigma_yy - sigma_zz)**2 + (sigma_zz - sigma_x)**2 ) / 2 + 3*((tau_xy))**2 + ((tau_yz))**2 + ((tau_xz))**2) )
    return