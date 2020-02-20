import stress_modules
import equilibrium

# TODO: Entire code when all stress modules and equilibrium are finished.


def some_function(length, parameter):
    # Dummy temp function
    return length+parameter


def displ(input_disp):
    # Adds up the displacements of all different inputs (stresses and equilibrium) on any axis required, outputs them
    # as a list of displacements
    n = len(input_disp)
    input_disp_a = some_function(n, input_disp[0])  # Equilibrium function
    input_disp_b = some_function(n, input_disp[1])  # Torsion function
    input_disp_c = some_function(n, input_disp[2])  # Bending function
    input_disp_d = some_function(n, input_disp[3])  # Stress function
    displace = []  # list of displacements
    for i in range(n):
        displace[i] = input_disp_a[i] + input_disp_b[i] + input_disp_c[i] + input_disp_d[i]
    return displace


def angle_disp(input_angle_disp):
    # Adds up all angle displacements from inputs (equilibrium, torsion, bending) and outputs them as a list
    n = len(input_angle_disp)
    input_angle_disp_a = some_function(n, input_angle_disp[0])
    input_angle_disp_b = some_function(n, input_angle_disp[1])
    input_angle_disp_c = some_function(n, input_angle_disp[2])
    tab = []
    for i in range(len(input_angle_disp)):
        tab[i] = input_angle_disp_a[i] + input_angle_disp_b[i] + input_angle_disp_c[i]
    return tab
