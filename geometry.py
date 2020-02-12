import numpy as np
from matplotlib import pyplot as plt
from math import sqrt, sin, cos, tan, acos


class Dataset:
    def __init__(self):
        # All the values are in SI units
        self.span = 1.691
        self.chord = 0.484
        self.hinge1 = 0.149
        self.hinge2 = 0.554
        self.hinge3 = 1.541
        self.height = 0.173
        self.skint = 0.0011
        self.spart = 0.0025
        self.stifft = 0.0012
        self.stiffh = 0.014
        self.stiffw = 0.018
        self.stiffn = 13
        self.xhinge1 = 0.149
        self.xhinge2 = 0.554
        self.xhinge3 = 1.541
        self.minx = - self.hinge2
        self.maxx = self.span - self.hinge2
        self.minz = -self.chord + self.height/2
        self.maxz = self.height/2
        self.miny = - self.height/2
        self.maxy = self.height/2

        self._radius = self.height*0.5
        self._a = sqrt(self._radius * self._radius + (self.chord-self._radius)*(self.chord-self._radius))
        self._circumference = 2*self._a + np.pi * self._radius

    def IzzSkinSpar(self):
        beta = acos(self.height*0.5/self._a)
        Izzskin1 = self.skint * self._a*self._a*self._a * sin(beta) * sin(beta)/6
        Izzskin2 = np.pi * self.skint * self.height*self.height*self.height / 16
        Izzspar = self.skint * self.height * self.height * self.height / 12
        return Izzskin1 + Izzskin2 + Izzspar

    def _stiffcoord(self, num):
        step = self._circumference/self.stiffn
        current = step*(num-1)
        if current < np.pi*0.25*self.height:
            angle = current / self._radius
            z = self._radius * cos(angle)
            y = - self._radius * sin(angle)
            return z, y
        elif current > self._circumference - np.pi*0.25*self.height:
            angle = (self._circumference - current) / self._radius
            z = self._radius * cos(angle)
            y = self._radius * sin(angle)
            return z, y
        elif current > np.pi*0.25*self.height + self._a:
            current = current - np.pi*0.25*self.height - self._a
            z = (self.chord - self._radius) * current/self._a - self.chord + self._radius
            y = self._radius * current/self._a
            return z, y
        elif current > np.pi*0.25*self.height:
            current -= np.pi*0.25*self.height
            z = - (self.chord - self._radius) * current/self._a
            y = - self._radius + self._radius * current/self._a
            return z, y
        else:
            raise ValueError("The aileron does not contain this number of stringers")

    def stiffLoc(self, n=None):
        if n == None:
            list = []
            for i in range(self.stiffn):
                list.append(self._stiffcoord(i))
            return list
        else:
            return self._stiffcoord(n)

    def visualinspection(self):
        plt.scatter(*zip(*self.stiffLoc()))
        plt.gca().invert_xaxis()
        plt.show()



