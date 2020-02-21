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
        # Comment because ???
        Iyy = 0
        zbar = self.centroid(2)
        if skin:
            beta = acos((self.chord - self.radius) / self.a)
            Iyy += self.skint * self.a * self.a * self.a * cos(beta) * cos(beta) * 2 / 3 + np.pi * self.skint * self.height * self.height * self.height / 16
            #Steiner term:
            Iyy += 2 * (self.skint * self.a) * ((self.chord - self.radius)*.5 - zbar )**2               # steiner terms for sloped part
            Iyy += (zbar + 2*self.radius / np.pi)**2 * np.pi * self.skint * self.radius                 # steiner terms for circular part
        if spar:
            Iyy += self.height * self.spart * self.spart * self.spart / 12
            # Steiner term:
            Iyy += self.spart * self.height * zbar**2                                                   # steiner term for spar
        if stiffener:
            Iyy += self._I_stiff(self.stiffLoc(), 0)
            # Steiner term:
            Iyy += sum([ (stiff[0]-self.radius - zbar)**2 * self.stiffener_area for stiff in self.stiffLoc()  ])    #steiner term for stiffeners
        return Iyy

    def centroid(self, axis=None):
        xbar = self.span * 0.5
        ybar = 0
        zbar = self.skint * (self.radius * self.radius * 2 - self.a * (self.chord - self.radius)) + (self.stiffener_area * sum(stiff[0] +0.0865 for stiff in self.stiffLoc()))
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

#    def centroid(self, axis=None):
#        xbar = self.span * 0.5
#        ybar = 0
#        zbar = self.skint * (self.radius ** 2 * 2 - self.a * (self.chord - self.radius)) + (self.stiffener_area * sum(stiff[0] +0.0865 for stiff in self.stiffLoc()))
#        zbar = zbar / (self.skint * (np.pi * self.radius * self.skint + 2 * self.a) + self.height * self.spart + self.stiffener_area * self.stiffn)
#        if axis == 0:
#            return xbar
#        if axis == 1:
#            return ybar
#        if axis == 2:
#            return zbar
#        if axis == None:
#            return xbar, ybar, zbar

    def shearcentre(self):
        xi = self.centroid(1)
        eta = 'centroid location in z'
        # We have the starting values of s for each section from 1 to 4, and we use the step value for each boom to
        # calculate their s position relative to s0 in each section. We can then use this to calculate each section
        # even if booms change in number, location or value. Currently not implemented and the implemented version
        # seems like it will be very fucking ugly, so that's to be fixed
        step = self._circumference/self.stiffn
        s_1 = 0
        s_2 = np.pi * self.radius * 0.5
        s_3 = s_2 + self.a
        s_4 = s_3 + self.a
        locations = self.stiffloc(shear=True)
        return xi, eta

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

    def stiffLoc(self, n=None, shear=False):
        # Returns a tuple of the given stiffener number.
        # If no arguments are given, a list of tuples containing all stiffener locations is returned.
        if n == None:
            list = []
            for i in range(self.stiffn):
                list.append(self._stiffcoord(i))
            return list
        elif shear:
            # Append to list elements in 4 divisions, signifying sections 1 through four in shear center calculation
            # Determine the s-value from
            list =[]
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
    """Class containing all the loads and tools for interpreting the applied loads.
    Can also calculate deflections within a discretized version of the aileron that has been idealized as a beam,
    using the aileron class.

    Initialization requires the number of sections that the aileron should be discretized in. Changing the number
    of sections should only be done before using any functions to calculate deflection, and only using the method
    self.gridnum(gridnum)."""
    def __init__(self, gridnum, filename="aerodata.csv", Nx=41, Nz=81, aileron=Aileron(), hinge1z=True,
                 hinge1d_y=0.00681, F_I_z=True, hinge2_fix=True, P=37900, hinge3z=True, hinge3d_y=0.0203, aero=True):
        self.filename = filename
        self.Nx = Nx
        self.Nz = Nz
        self.a = aileron
        self.aerogrid = self.aero_points(self.Nx, self.Nz, self.a)
        self.res_locs, self.res_forces = self.get_aero_resultants(self.filename, self.aerogrid)
        self.q = partial(self._getq, self.aerogrid, self.res_forces, self.a)
        self._grid = self.make_cuts_sections(self.a, gridnum)
        self.gridsize = self._grid.size
        self.dx_list = np.asarray([self._grid[i + 1] - self._grid[i] for i in range(self._grid.size - 1)])
        self.hinge1_idx = np.argmin(np.abs(self._grid - self.a.hinge1))
        self.hinge1_val = self._grid[self.hinge1_idx]
        self.hinge2_idx = np.argmin(np.abs(self._grid - self.a.hinge2))
        self.hinge2_val = self._grid[self.hinge2_idx]
        self.hinge3_idx = np.argmin(np.abs(self._grid - self.a.hinge3))
        self.hinge3_val = self._grid[self.hinge3_idx]

        self._geo_error = np.min(self.dx_list)*0.1

        self.hinge1z = hinge1z
        self.hinge1d_y = hinge1d_y
        self.hinge3z = hinge3z
        self.hinge3d_y = hinge3d_y
        self.P = P
        self.jammed = F_I_z
        self.hinge2_fix = hinge2_fix

        if self.hinge1d_y != 0 and self.hinge3d_y != 0:
            self.bending=True
        else:
            self.bending=False

    def renew_grid(self, gridnum):
        self._grid = self.make_cuts_sections(self.a, gridnum)
        self.dx_list = [self._grid[i + 1] - self._grid[i] for i in range(self._grid.size - 1)]
        return None

    def gridnum(self, gridnum):
        self.renew_grid(gridnum)
        return None

    def grid(self, i):
        return self._grid[i]

    def hinge_idx(self, num):
        if num == 1:
            return self.hinge1_idx
        elif num == 2:
            return self.hinge2_idx
        elif num == 3:
            return self.hinge3_idx
        else:
            raise ValueError("This hinge does not exist.")

    def hinge_val(self, num):
        if num == 1:
            return self.hinge1_val
        elif num == 2:
            return self.hinge2_val
        elif num == 3:
            return self.hinge3_val
        else:
            raise ValueError("This hinge does not exist.")

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


    def _getq(self, coords, forces, a, x):
        xlocs = np.unique(coords[:, 0])
        if x < xlocs[0]:
            return forces[0]
        elif x > xlocs[-1]:
            return forces[-1]
        return interpolate(xlocs, forces, x)


    def make_cuts_sections(self, a, num_sections_x):
        new_list = [0]
        step = a.span/num_sections_x
        for i in range(1,num_sections_x):
            new_list.append(step*i)
        new_list.append(a.span)
        return np.asarray(new_list)

    def distr_defl_func(self, q1, q2, L, I, E=71000000000):
        # deflection return positive in positive direction of q
        if q1 != q2:
            return L*L*L*L/(E*I)*(q2 /8 + (q1-q2)*3/72)
        else:
            return L*L*L*L*q2/(E*I*8)


    def distr_angle_func(self, q1, q2, L, I, E=71000000000):
        if q1 != q2:
            return L*L*L/(E*I)*(q2/6 + (q1-q2)/18)
        else:
            return L*L*L*q2/(6*E*I)

    def pointload_defl(self, P, L, I, E=71000000000):
        return P*L*L*L/(3*E*I)

    def pointload_angle(self, P, L, I, E=71000000000):
        return P*L*L/(2*E*I)

    def forward_defl_point(self, i):
        if abs(self._grid[i]-self.hinge1_val) < self._geo_error:
            return self.pointload_defl(1, self.dx_list[i], self.a.Izz())
        else:
            return 0

    def backward_defl_point(self, i):
        if abs(self._grid[i]-self.hinge3_val) < self._geo_error:
            return self.pointload_defl(1, self.dx_list[i-1], self.a.Izz())
        else:
            return 0

    def forward_angle_point(self, i):
        if abs(self._grid[i]-self.hinge1_val) < self._geo_error:
            return self.pointload_angle(1, self.dx_list[i], self.a.Izz())
        else:
            return 0

    def backward_angle_point(self, i):
        if abs(self._grid[i]-self.hinge3_val) < self._geo_error:
            return self.pointload_angle(1, self.dx_list[i-1], self.a.Izz())
        else:
            return 0

    def forward_defl_distr(self, i):
        return self.distr_defl_func(-self.q(self._grid[i]), -self.q(self._grid[i + 1]), self.dx_list[i], self.a.Izz())

    def backward_defl_distr(self, i):
        return self.distr_defl_func(-self.q(self._grid[i - 1]), -self.q(self._grid[i]), self.dx_list[i - 1], self.a.Izz())

    def forward_angle_distr(self, i):
        return self.distr_angle_func(-self.q(self._grid[i]), -self.q(self._grid[i + 1]), self.dx_list[i], self.a.Izz())

    def backward_angle_distr(self, i):
        return self.distr_angle_func(-self.q(self._grid[i - 1]), -self.q(self._grid[i]), self.dx_list[i - 1], self.a.Izz())




if __name__ == "__main__":
    print("Hello world")
