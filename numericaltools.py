import numpy as np

#Integration function, with func being the (mathematical) function to integrate, start and stop being beginning and end
#of the integration respectively, and number_of_points the subdivisions.
def integrate(func, start, stop, number_of_points):
    _integration_temp = 0
    _diff = (start+stop)/number_of_points
    for N in range(0, number_of_points):
        _point_a = start+N*_diff
        _point_b = _point_a+_diff
        _integration_temp += _diff/2*(func(_point_a)+func(_point_b))
    integration = _integration_temp
    return integration


# function which calculates the value of each section in the linear spline,
# also can be referenced as a spline item.
def splineitem(fi, f1, xi, x1, x):
    #
    splinei = fi + ((f1-fi)/(x1-xi))*(x-xi)
    return splinei


def spline(x):
    #   Write code which for every range x-xi calls upon the splineitem function
    #   to build a linear spline.
    sp = 1 * x
    return sp
