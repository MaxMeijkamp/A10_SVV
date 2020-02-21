import numpy as np
from matplotlib import pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from InputClasses import Aileron
from math import sin, cos
from numericaltools import integrate, interpolate, cont_spline
from functools import *


def get_theta(i, N):
    return (i-1)*np.pi/N


def getzcoord(i, aileron, Nz=81):
    return -0.5 * (aileron.chord*0.5 * (1-cos(get_theta(i, Nz))) + aileron.chord*0.5*(1-cos(get_theta(i+1, Nz))))


def getxcoord(i, aileron, Nx=41):
    return 0.5 * (aileron.span*0.5 * (1-cos(get_theta(i, Nx))) + aileron.span*0.5*(1-cos(get_theta(i+1, Nx))))


def aero_points(Nx, Nz, a):
    coordlist = []
    for x in range(1, Nx+1):
        for z in range(1, Nz+1):
            coordlist.append([getxcoord(x, a), getzcoord(z, a)])
    return np.asarray(coordlist)


def make_sections(Nx, Nz, a):
    coords = np.unique(aero_points(Nx, Nz, a)[:, 0])
    new_list = [0]
    for i in range(coords.size-1):
        new_list.append(0.5*(coords[i] + coords[i+1]))
    new_list.append(a.span)
    return np.asarray(new_list)


def get_aero_resultants(file, coords):
    data = np.genfromtxt(file, delimiter=",")
    data[0,0] = 0.034398 # Reader gives the first value as nan, this is to fix that problem
    coords = np.unique(coords[:,1])
    res_forces = []
    res_locations = []
    data = np.flip(data, axis=0)
    for i in range(data.shape[1]):
        Q = 0
        res_forces.append(integrate(cont_spline(coords, data[:,i]), np.min(coords), np.max(coords), 100))
        for j in range(data.shape[0]):
            Q += data[j,i]*coords[j]
        res_locations.append(Q/np.sum(data[:,i]))
    return res_locations, res_forces

def getq(coords, forces, a, x):
    if x == 0:
        return forces[0]
    elif x == a.span:
        return forces[-1]
    xlocs = np.unique(coords[:,0])
    return interpolate(xlocs, forces, x)


if __name__ == "__main__":
    a = Aileron()
    aerogrid = aero_points(41, 81, a)
    locs, forces = get_aero_resultants("aerodata.csv", aerogrid)
    print(np.unique(aerogrid[:,0]))
    print(forces)
    print(getq(aerogrid, forces, a, 0.5))

    # plt.scatter(aerogrid[:,0], aerogrid[:,1], s=10)
    # plt.scatter(np.unique(aerogrid[:,0]), locs, s=50)
    # plt.gca().invert_xaxis()
    # plt.plot([0, a.span], [-a.chord, -a.chord])
    # plt.plot([0, a.span], [0, 0])
    # plt.plot([0, 0], [-a.chord, 0])
    # plt.plot([a.span, a.span], [-a.chord, 0])
    # plt.scatter(a.hinge2, -a._radius, s=100)
    # plt.scatter(a.hinge1, -a._radius, s=100)
    # plt.scatter(a.hinge3, -a._radius, s=100)
    # plt.show()
    #
    #
    #
    #
    #

