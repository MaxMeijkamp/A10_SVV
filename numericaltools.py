import numpy as np
from functools import *


def integrate(func, start, stop, number_of_points):
    # Integration function, with func being the (mathematical) function to integrate,
    # start and stop being beginning and end of the integration respectively,
    # and number_of_points is the subdivisions. (Using Numpy)
    start, stop, number_of_points = float(start), float(stop), float(number_of_points)
    width = (stop - start) / number_of_points
    points = np.linspace(start, stop - width, number_of_points)
    next_points = np.linspace(start + width, stop, number_of_points)
    values, next_values = func(points), func(next_points)
    integration = (np.sum((next_values - values) * .5) + np.sum(values)) * width
    return integration


def cont_spline(x_discrete, f_discrete):
    return np.vectorize(partial(interpolate, x_discrete, f_discrete))


def spline(x, f, n):
    # Spline function,
    if len(x) != len(f):
        raise ValueError("The lists are not of the same shape")
    sp_start = []
    sp_slope = []
    for i in range(0, n-1):
        splinestart = f[i]
        splineslope = ((f[i+1] - f[i])/(x[i+1] - x[i]))
        sp_start.append(splinestart)
        sp_slope.append(splineslope)
    # sp = np.vstack((np.array(sp_start), np.array(sp_slope)))
    return sp_start, sp_slope  # sp.T


def interpolate(x, f, x_target):
    # Interpolation function which gives f_target for a given x_target and x-f spline
    n = len(x)
    if x[0] > x_target or x_target > x[n-1]:
        raise ValueError("The target location is out of bounds")
    sp_start, sp_slope = spline(x, f, n)
    left_i = n-2
    for i in range(n):
        if x[i] > x_target:
            left_i = i-1
            break
    f_target = sp_slope[left_i] * (x_target - x[left_i]) + sp_start[left_i]
    return f_target
