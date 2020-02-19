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

#Integration function, with func being the (mathematical) function to integrate, start and stop being beginning and end
#of the integration respectively, and number_of_points is the subdivisions. (Using Numpy)
def integrateP(func, start, stop, number_of_points):
    start, stop, number_of_points = float(start), float(stop), float(number_of_points)
    width = (stop - start) / (number_of_points)
    points = np.linspace(start, stop -  width, number_of_points)
    next_points = np.linspace(start +  width, stop, number_of_points)
    values, next_values = func(points), func(next_points)
    integration = (np.sum((next_values - values) * .5) + np.sum(values)) * width
    return integration
