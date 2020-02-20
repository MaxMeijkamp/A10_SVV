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

def findShearCenter(dataset):
    Izz = dataset.Izz()
    SC_y = dataset.centroid[1]
    Vy = 1 #force applied in shear centre
    #cuts at front of aileron and middle in spar
    #y,z axis on middle of spar
    r = dataset._radius
    a = dataset.a
    skint = dataset.skint
    height = dataset.maxy
    stiffenerA = dataset._stiffener_area
    #funcs will return y as function of s
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
    booms_s1 = [boom for boom in dataset.stiffLoc() if boom[1] >= 0 and boom[0] > 0]
    booms_s2 = [boom for boom in dataset.stiffLoc() if boom[1] >= 0 and boom[0] < 0]
    booms_s3 = [boom for boom in dataset.stiffLoc() if boom[1] < 0 and boom[0] < 0]
    booms_s4 = [boom for boom in dataset.stiffLoc() if boom[1] < 0 and boom[0] > 0]
    booms_s5 = [boom for boom in dataset.stiffLoc() if boom[1] < 0 and boom[0] == 0]
    booms_s6 = [boom for boom in dataset.stiffLoc() if boom[1] >= 0 and boom[0] == 0]
    all_booms_divided = [booms_s1, booms_s2, booms_s3, booms_s4, booms_s4, booms_s5, booms_s6]

    for boomlist in all_booms_divided:
        for boom in boomlist:
            boom = boom[1] * stiffenerA

    dqb1 = -1. / Izz * (numericaltools.integrate(func_s1, 0, dist_s1, 100000)  +  sum(all_booms_divided[0]) )
    dqb2 = -1. / Izz * (numericaltools.integrate(func_s2, 0, dist_s2, 100000)  +  sum(all_booms_divided[1]) )
    dqb3 = -1. / Izz * (numericaltools.integrate(func_s3, 0, dist_s3, 100000)  +  sum(all_booms_divided[2]) )
    dqb4 = -1. / Izz * (numericaltools.integrate(func_s4, 0, dist_s4, 100000)  +  sum(all_booms_divided[3]) )
    dqb5 = -1. / Izz * (numericaltools.integrate(func_s5, 0, dist_s5, 100000)  +  sum(all_booms_divided[4]) )
    dqb6 = -1. / Izz * (numericaltools.integrate(func_s6, 0, dist_s6, 100000)  +  sum(all_booms_divided[5]) )

    # now find the qs,0 by setting the twist to 0: ie 1/2A * integral(q ds / (G t)) = 0
    # G = constant, A is constant per section. t = constant and q = sum of (qb times the length - qs, 0)
    # therefor the qs,0 should be opposite of the sum of qb
    # clockwise positive

    #for left part:

    #round part:
    dqb1 +=5


    return 0

