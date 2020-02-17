import numpy as np


# Integration function, with func being the (mathematical) function to integrate, start and stop being beginning and end
# of the integration respectively, and number_of_points the subdivisions.
def integrate(func, start, stop, number_of_points):
    integration_temp = 0
    diff = (start+stop) / number_of_points
    for N in range(0, number_of_points):
        point_a = start + N*diff
        point_b = point_a + diff
        integration_temp += diff / 2 * (func(point_a) + func(point_b))
    integration = integration_temp
    return integration


# Spline function,
def spline(x, f, n):
    sp_start = []
    sp_slope = []
    for i in range(0,n):
        splinestart = f[i] + ((f[i+1] - f[i])/(x[i+1] - x[i]))*x[i]
        splineslope = ((f[i+1] - f[i])/(x[i+1] - x[i]))
        sp_start.append(splinestart)
        sp_slope.append(splineslope)
    sp = np.empty(shape=[n, 2])
    return sp
