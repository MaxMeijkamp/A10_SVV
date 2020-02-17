import geometry

def bending(y, z, Ixx, Izz, Iyy, M_z, M_y, centroid):
    #returns sigma_x for location y, z on the aileron. Inputs: the point x, y;
    #Cont'd: the MOI Izz, Iyy; the relevant moments M_z, M_y (tension in + plane positive), centroid: tuple with xyz centroid of the aileron
    return M_y *  (z - centroid[2]) / Iyy + M_z * (y - centroid[1]) / Izz