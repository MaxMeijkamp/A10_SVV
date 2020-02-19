import stress_modules
import equilibrium


def some_function(n, l):
    # Dummy temp function
    return n+l


def displ(input):
    # Adds up the displacements of all different inputs (stresses and equilibrium) on any axis required, outputs them
    # as a list of displacements
    n = len(input)
    input_disp_a = some_function(n, a)  # Equilibrium function
    input_disp_b = some_function(n, b)  # Torsion function
    input_disp_c = some_function(n, c)  # Bending function
    input_disp_d = some_function(n, c)  # Stress function
    displace = []  # list of displacements
    for i in range(n):
        displace[i] = input_disp_a[i] + input_disp_b[i] + input_disp_c[i] + input_disp_d[i]
    return displace


def angle_disp(input_angle_disp):
    # Adds up all angle displacements from inputs (equilibrium, torsion, bending) and outputs them as a list
    n = len(input)
    input_angle_disp_a = some_function(n, a)
    input_angle_disp_b = some_function(n, b)
    input_angle_disp_c = some_function(n, c)
    tab = []
    for i in range(len(input_angle_disp)):
        tab[i] = input_angle_disp_a[i] + input_angle_disp_b[i] + input_angle_disp_c[i]
    return tab
