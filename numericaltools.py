import numpy as np


def integrate(func, start, stop, number_of_points):
    integration_temp = 0
    diff = (start+stop) / number_of_points
    for N in range(0, number_of_points):
        point_a = start + N*diff
        point_b = point_a + diff
        integration_temp += diff / 2 * (func(point_a) + func(point_b))
    integration = integration_temp
    return integration


def spline(x, f):
    # Spline function,
    if len(x) != len(f):
        raise ValueError("The lists are not of the same shape")
    n = len(x)
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
    sp_start, sp_slope = spline(x, f)
    for i in range(0, len(x)):
        if x[i] >= x_target:
            left_i = i-1
            break
    f_target = sp_slope[left_i] * (x_target - x[left_i]) + sp_start[left_i]
    return f_target


def integrateP(func, start, stop, number_of_points):
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


if __name__ == "__main__":
    x = np.array()
    f = np.array()
    interp = interpolate([0, 0.1, 0.2, 0.5, 0.65, 0.8, 0.9, 1.0], [1, 4, 2, 3, 7, -4, 0.5, 0], 0.15)
    print(interp)
