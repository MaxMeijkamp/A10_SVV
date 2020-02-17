import numpy as np
from matplotlib import pyplot as plt
from math import sqrt, sin, cos, tan, acos


class Dataset:
    def __init__(self, span=1.691, chord=0.484, hinge1=0.149, hinge2=0.554, hinge3=1.541, height=0.173, skint=0.0011, spart=0.0025, stifft=0.0012, stiffh=0.014, stiffw=0.018, stiffn=13):
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

        # Protected variables
        self._radius = self.height*0.5
        self._a = sqrt(self._radius * self._radius + (self.chord-self._radius)*(self.chord-self._radius))
        self._circumference = 2*self._a + np.pi * self._radius
        self._stiffener_area = self.stifft*(self.stiffh + self.stiffw)

    def Izz(self, skin=True, spar=True, stiffener=True):
        # Calculates Izz of a cross-section. Also able to calculate only parts of Izz based on arguments given
        Izz = 0
        if skin:
            beta = acos((self.chord - self._radius) / self._a)
            Izz += self.skint * self._a*self._a*self._a * sin(beta) * sin(beta) * 2/3 + np.pi * self.skint * self.height*self.height*self.height / 16
        if spar:
            Izz += self.skint * self.height * self.height * self.height / 12
        if stiffener:
            Izz += self._Izzstiff(self.stiffLoc())
        return Izz

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
        # Returns a tuple of the given stiffener number. If no arguments are given, a list of tuples containing all stiffener locations is returned.
        if n == None:
            list = []
            for i in range(self.stiffn):
                list.append(self._stiffcoord(i))
            return list
        else:
            return self._stiffcoord(n)

    def _Izzstiff(self, stiffener_list):
        Izz = 0
        for stiff in stiffener_list:
            d_2 = stiff[1]*stiff[1]
            Izz += d_2*self._stiffener_area
        return Izz

    def visualinspection(self):
        # Used as a visual inspection of the cross-section.
        plt.scatter(*zip(*self.stiffLoc()))
        plt.gca().invert_xaxis()
        y1 = np.linspace(0,np.pi,1000)
        z1 = np.sin(y1)*self._radius
        y1 = np.cos(y1)*self._radius
        y2 = self._radius * np.linspace(-1,1,1000)
        z2 = np.linspace(0,-self.chord+self._radius, 500)
        z2 = np.append(z2, -z2-self.chord+self._radius)
        z3 = np.zeros(100)
        y3 = np.linspace(-self._radius, self._radius, 100)
        plt.plot(z1, y1)
        plt.plot(z2, y2)
        plt.plot(z3, y3)
        plt.axis('equal')
        plt.xlabel("$z \: [m]$")
        plt.ylabel("$y \: [m]$")
        plt.show()



