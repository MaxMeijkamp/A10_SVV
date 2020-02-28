import numpy as np
from math import sin, cos
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

def calc_shearstress(a, x, x_pos, Sy, Sz, My, Mz, T):
    idx = np.where(x==x_pos)[0]
    cutoff3 = a._circumference - np.pi * 0.25 * a.height
    cutoff2 = np.pi * 0.25 * a.height + a.a
    cutoff1 = np.pi * 0.25 * a.height
    points = np.linspace(0, a._circumference, 1000)
    points3 = points[np.where(points>cutoff3)]
    points2 = points[np.where(points>cutoff2)]
    points2 = points2[np.where(points2<cutoff3)]
    points1 = points[np.where(points>cutoff1)]
    points1 = points1[np.where(points1<cutoff2)]
    points0 = points[np.where(points>0)]
    points0 = points0[np.where(points0<cutoff1)]

    point_list = np.array([])
    for point in points3:
        angle = (a._circumference - point) / a.radius
        z = a.radius * cos(angle)
        y = a.radius * sin(angle)
    if current < np.pi * 0.25 * self.height:
        angle = current / self.radius
        z = self.radius * cos(angle)
        y = - self.radius * sin(angle)
        return z, y
    elif current > self._circumference - np.pi * 0.25 * self.height:
        angle = (self._circumference - current) / self.radius
        z = self.radius * cos(angle)
        y = self.radius * sin(angle)
        return z, y
    elif current > np.pi * 0.25 * self.height + self.a:
        current = current - np.pi * 0.25 * self.height - self.a
        z = (self.chord - self.radius) * current / self.a - self.chord + self.radius
        y = self.radius * current / self.a
        return z, y
    elif current > np.pi * 0.25 * self.height:
        current -= np.pi * 0.25 * self.height
        z = - (self.chord - self.radius) * current / self.a
        y = - self.radius + self.radius * current / self.a
        return z, y
    else:
        raise ValueError("The aileron does not contain this number of stiffeners")

    sigma_11 = bending([My[idx], Mz[idx]], point_vect, a)

def bending(M_vect, point_vect, a):
    # returns sigma_x for location y, z on the aileron. Inputs: Internal moment vector (y, z); position of the point to
    # calculate the stresses (y, z); an instance of Aileron class.
    # the relevant moments M_z, M_y (tension in + plane positive), centroid: tuple with xyz centroid of the aileron
    return - M_vect[1] * (point_vect[0] - a.centroid(axis=1)) / a.Izz() + M_vect[0] * (point_vect[1] - a.centroid(2)) / a.Iyy()

# For bending, not only sigma_x is needed, but also the displacements and angular displacements
# at all sections caused by the bending.
