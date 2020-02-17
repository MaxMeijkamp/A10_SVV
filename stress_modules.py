import numpy as np

def von_mises()
    #Takes arg. simga_xx, sigma_yy, sigma_zz, tau_xy, tau_yz, tau_xz
    sigma_vm = np.sqrt( ( (sigma_xx - sigma_yy)**2 + (sigma_yy - sigma_zz)**2 + (sigma_zz - sigma_x)**2 ) / 2 + 3*((tau_xy))**2 + ((tau_yz))**2 + ((tau_xz))**2) )
    return