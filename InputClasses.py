import numpy as np
from matplotlib import pyplot as plt
from math import sqrt, sin, cos, acos
from numericaltools import *


class Aileron:
    def __init__(self, span=1.691, chord=0.484, hinge1=0.149, hinge2=0.554, hinge3=1.541, height=0.173, skint=0.0011,
                 spart=0.0025, stifft=0.0012, stiffh=0.014, stiffw=0.018, stiffn=13):
        # All given parameters and useful parameters. All the values are in SI units
        self.span = span
        self.chord = chord
        self.hinge1 = hinge1
        self.hinge2 = hinge2
        self.hinge3 = hinge3
        self.height = height
        self.skint = skint
        self.spart = spart
        self.stifft = stifft
        self.stiffh = stiffh
        self.stiffw = stiffw
        self.stiffn = stiffn
        self.xhinge1 = self.hinge1-self.hinge2
        self.xhinge2 = 0.0
        self.xhinge3 = self.hinge3-self.hinge2
        self.minx = - self.hinge2
        self.maxx = self.span - self.hinge2
        self.minz = -self.chord + self.height/2
        self.maxz = self.height/2
        self.miny = - self.height/2
        self.maxy = self.height/2
        self.stiffener_area = self.stifft * (self.stiffh + self.stiffw)
        self.radius = self.height * 0.5
        self.a = sqrt(self.radius * self.radius + (self.chord - self.radius) * (self.chord - self.radius))

        # Protected variables
        self._circumference = 2 * self.a + np.pi * self.radius

    def Izz(self, skin=True, spar=True, stiffener=True):
        # Calculates Izz of a cross-section. Also able to calculate only parts of Izz based on arguments given
        Izz = 0
        if skin:
            beta = acos((self.chord - self.radius) / self.a)
            Izz += self.skint * self.a * self.a * self.a * sin(beta) * sin(beta) * 2 / 3 + np.pi * self.skint * self.height * self.height * self.height / 16
        if spar:
            Izz += self.spart * self.height * self.height * self.height / 12
        if stiffener:
            Izz += self._I_stiff(self.stiffLoc(), 1)
        return Izz

    def Iyy(self, skin=True, spar=True, stiffener=True):
        # Calculates Iyy of a cross-section. Also able to calculate only parts of Iyy based on arguments given
        Iyy = 0
        if skin:
            beta = acos((self.chord - self.radius) / self.a)
            Iyy += self.skint * self.a * self.a * self.a * cos(beta) * cos(beta) / 12 + np.pi * self.skint * self.height * self.height * self.height / 16
        if spar:
            Iyy += self.height * self.spart * self.spart * self.spart / 12
        if stiffener:
            Iyy += self._I_stiff(self.stiffLoc(), 0)
        return Iyy

    def centroid(self, axis=None):
        xbar = self.span * 0.5
        ybar = 0
        zbar = self.skint * (self.radius * self.radius * 2 - self.a * (self.chord - self.radius)) + sum([self.stiffener_area * stiff[0] for stiff in self.stiffLoc()])
        zbar = zbar / (self.skint * self._circumference + self.spart * self.height + self.stiffener_area * self.stiffn)
        if axis == 0:
            return xbar
        if axis == 1:
            return ybar
        if axis == 2:
            return zbar
        if axis == None:
            return xbar, ybar, zbar
        else:
            raise ValueError("The axis is invalid")

    def shearcentre(self, axis=None):
        # TODO: implement
        return self.centroid(axis)

    def _stiffcoord(self, num):
        step = self._circumference/self.stiffn
        current = step*(num-1)
        if current < np.pi*0.25*self.height:
            angle = current / self.radius
            z = self.radius * cos(angle)
            y = - self.radius * sin(angle)
            return z, y
        elif current > self._circumference - np.pi*0.25*self.height:
            angle = (self._circumference - current) / self.radius
            z = self.radius * cos(angle)
            y = self.radius * sin(angle)
            return z, y
        elif current > np.pi*0.25*self.height + self.a:
            current = current - np.pi*0.25*self.height - self.a
            z = (self.chord - self.radius) * current / self.a - self.chord + self.radius
            y = self.radius * current / self.a
            return z, y
        elif current > np.pi*0.25*self.height:
            current -= np.pi*0.25*self.height
            z = - (self.chord - self.radius) * current / self.a
            y = - self.radius + self.radius * current / self.a
            return z, y
        else:
            raise ValueError("The aileron does not contain this number of stringers")

    def stiffLoc(self, n=None):
        # Returns a tuple of the given stiffener number.
        # If no arguments are given, a list of tuples containing all stiffener locations is returned.
        if n == None:
            list = []
            for i in range(self.stiffn):
                list.append(self._stiffcoord(i))
            return list
        else:
            return self._stiffcoord(n)

    def _I_stiff(self, stiffener_list, axis):
        # Calculate moment of Inertia, including Steiner terms
        i_stiff = 0
        for stiff in stiffener_list:
            d_2 = stiff[axis]*stiff[axis]
            i_stiff += d_2*self.stiffener_area
        return i_stiff

    def visualinspection(self):
        # Used as a visual inspection of the cross-section.
        plt.scatter(*zip(*self.stiffLoc()))
        plt.gca().invert_xaxis()
        y1 = np.linspace(0, np.pi, 1000)
        z1 = np.sin(y1)*self.radius
        y1 = np.cos(y1)*self.radius
        y2 = self.radius * np.linspace(-1, 1, 1000)
        z2 = np.linspace(0, -self.chord + self.radius, 500)
        z2 = np.append(z2, -z2 - self.chord + self.radius)
        z3 = np.zeros(100)
        y3 = np.linspace(-self.radius, self.radius, 100)
        plt.plot(z1, y1)
        plt.plot(z2, y2)
        plt.plot(z3, y3)
        plt.scatter(self.centroid(2), self.centroid(1), s=100)
        plt.axis('equal')
        plt.xlabel("$z \: [m]$")
        plt.ylabel("$y \: [m]$")
        plt.show()


class AppliedLoads:
    def __init__(self, filename="aerodata.csv", Nx=41, Nz=81, aileron=Aileron()):
        self.filename = filename
        self.Nx = Nx
        self.Nz = Nz
        self.a = aileron
        self.aerogrid = self.aero_points(self.Nx, self.Nz, self.a)
        self.res_locs, self.res_forces = self.get_aero_resultants(self.filename, self.aerogrid)

    def _getzcoord(self, i, aileron, Nz=81):
        return -0.5 * (aileron.chord*0.5 * (1-cos(self._get_theta(i, Nz))) + aileron.chord*0.5*(1-cos(self._get_theta(i+1, Nz))))

    def _getxcoord(self, i, aileron, Nx=41):
        return 0.5 * (aileron.span*0.5 * (1-cos(self._get_theta(i, Nx))) + aileron.span*0.5*(1-cos(self._get_theta(i+1, Nx))))

    def _get_theta(self, i, N):
        return (i - 1) * np.pi / N

    def aero_points(self, Nx, Nz, a):
        coordlist = []
        for x in range(1, Nx+1):
            for z in range(1, Nz+1):
                coordlist.append([self._getxcoord(x, a), self._getzcoord(z, a)])
        return np.asarray(coordlist)

    def get_aero_resultants(self, file, coords):
        data = np.genfromtxt(file, delimiter=",")
        data[0, 0] = 0.034398  # Reader gives the first value as nan, this is to fix that problem
        coords = np.unique(coords[:, 1])
        res_forces = []
        res_locations = []
        data = np.flip(data, axis=0)
        for i in range(data.shape[1]):
            Q = 0
            res_forces.append(integrate(cont_spline(coords, data[:, i]), np.min(coords), np.max(coords), 100))
            for j in range(data.shape[0]):
                Q += data[j, i] * coords[j]
            res_locations.append(Q / np.sum(data[:, i]))
        return res_locations, res_forces



if __name__ == "__main__":
    print("Hello world")
