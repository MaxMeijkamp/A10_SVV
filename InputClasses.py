import numpy as np
from matplotlib import pyplot as plt
from math import sqrt, sin, cos, acos
from numericaltools import *

class Aileron:
    def __init__(self, span=1.691, chord=0.484, hinge1=0.149, hinge2=0.554, hinge3=1.541, height=0.173, skint=0.0011,
                 spart=0.0025, stifft=0.0012, stiffh=0.014, stiffw=0.018, stiffn=13, actdist=0.272):
        # All given parameters and useful parameters. All the values are in SI units
        self.span = span
        self.chord = chord
        self.hinge1 = hinge1
        self.hinge2 = hinge2
        self.hinge3 = hinge3
        self.actdist = actdist
        self.act1 = self.hinge2-0.5*self.actdist
        self.act2 = self.hinge2+0.5*self.actdist
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
        self.xact1 = self.act1-self.hinge2
        self.xact2 = self.act2-self.hinge2
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
            Izz += self.skint * self.a * self.a * self.a * sin(beta) * sin(beta) * 2 / 3 + np.pi / 2 * self.radius * self.radius * self.radius * self.skint
        if spar:
            Izz += self.spart * self.height * self.height * self.height / 12
        if stiffener:
            Izz += self._I_stiff(self.stiffLoc(), 1)
        return Izz

    def Iyy(self, skin=True, spar=True, stiffener=True):
        # Calculates Iyy of a cross-section. Also able to calculate only parts of Iyy based on arguments given
        Iyy = 0
        zbar = self.centroid(2)
        if skin:
            beta = acos((self.chord - self.radius) / self.a)
            Iyy += (np.pi/2 - 4/np.pi) * self.skint * self.radius * self.radius * self.radius
            Iyy += 2 * (self.skint * self.a * self.a * self.a * cos(beta) * cos(beta) / 12)
            Iyy += 2 * self.skint * self.a * ((-self.chord + self.radius) * .5 - zbar) * ((-self.chord + self.radius) * .5 - zbar)
            Iyy += (-zbar + 2*self.radius / np.pi) * (-zbar + 2*self.radius / np.pi) * np.pi * self.skint * self.radius
        if spar:
            Iyy += self.spart * self.height * zbar* zbar
        if stiffener:
            Iyy += self._I_stiff(self.stiffLoc(), 0)
            return Iyy

    def centroid(self, axis=None):
        xbar = self.span * 0.5
        ybar = 0
        zbar = self.skint * (self.radius * self.radius * 2 - self.a * (self.chord - self.radius))
        zbar += (self.stiffener_area * sum(stiff[0] for stiff in self.stiffLoc()))
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

    def shearcentre(self, verif = True):
        if verif:
            return -0.00535594953
        r = self.radius
        a = self.a
        Izz = self.Izz()
        skint = self.skint

        dist_s1 = dist_s4 = .5 * np.pi * r
        dist_s2 = dist_s3 = a
        dist_s5 = dist_s6 = r

        def func_s(s):
            if s <= dist_s1:
                return (r * np.sin(s / r)) * skint
            elif dist_s1 < s <= dist_s1 + dist_s2:
                return (r - (s - dist_s1) * r / a) * skint
            elif dist_s2 + dist_s1 < s <= dist_s3 + dist_s2 + dist_s1:
                return (0. - (s - dist_s1 - dist_s2) * r / a) * skint
            elif dist_s3 + dist_s2 + dist_s1 < s <= dist_s4 + dist_s3 + dist_s2 +dist_s1:
                return (r * -np.cos((s - .5 * r * np.pi - a * 2) / r)) * skint

        def func_s56(s):  # from s = 0 to s = height / 2
            return s * self.spart

        integrateV = np.vectorize(integrate)
        func_s = np.vectorize(func_s)
        num_steps = 10000

        s_list = np.zeros((7, num_steps)) # array with s value, integral shear, booms
        s_list[0] = np.linspace(0, self._circumference, num_steps)
        idx_s = [
                0, int((np.abs(s_list[0] - dist_s1)).argmin()),
                int(np.abs(s_list[0] - dist_s1 - dist_s2).argmin()),
                int(np.abs(s_list[0] - dist_s1 - dist_s2 - dist_s3).argmin()),
                int(np.abs(s_list[0] - dist_s1 - dist_s2 - dist_s3 - dist_s4).argmin())
                ]

        ss = np.arange(-r, r, self._circumference / num_steps)
        s56_list = np.zeros((7, ss.size))

        s56_list[0] = ss
        idx_s56 = [int(ss.size/2), int(ss.size-1)]
        s56_list[1] = - 1. / Izz * (integrateV(func_s56, 0, s56_list[0], 500))
        '''
        # Base shear flow of skin:
        s_list[1, 0: idx_s[1]] = -1. / Izz * (integrateV(func_s, 0, s_list[0, 0 : idx_s[1]], 500))
        s_list[1, idx_s[1]: idx_s[2]] = -1. / Izz * (integrateV(func_s, s_list[0, idx_s[1]], s_list[0, idx_s[1]: idx_s[2]], 500)) + s_list[1, idx_s[1]-1]
        s_list[1, idx_s[2]: idx_s[3]] = -1. / Izz * (integrateV(func_s, s_list[0, idx_s[2]], s_list[0, idx_s[2]: idx_s[3]], 500)) + s_list[1, idx_s[2]-1]
        s_list[1, idx_s[3]: idx_s[4]] = -1. / Izz * (integrateV(func_s, s_list[0, idx_s[3]], s_list[0, idx_s[3]: idx_s[4]], 500)) + s_list[1, idx_s[3]-1]
        s_list[1, idx_s[4]] = -1. / Izz * (integrateV(func_s, s_list[0, idx_s[3]], s_list[0, idx_s[4]], 500)) + s_list[1, idx_s[3]-1]
        '''

        s_list[1] = -1. / Izz * (integrateV(func_s, 0, s_list[0], 500))

        #Base shear flow of booms:
        stepsize = self._circumference/self.stiffn

        for i in range(self.stiffn):
            boomdist = stepsize * i
            idx = (np.abs(s_list[0] - boomdist)).argmin()
            s_list[2, idx : -1] +=  -1. / Izz * self.stiffener_area * func_s(s_list[0, idx]) / self.skint

        # Combined in row 3:
        s_list[3] = s_list[2] + s_list[1]

        stepcorrection = (num_steps - 1) / num_steps
        # Find correction shear flow:
        stepwidth = self._circumference / num_steps


        splineshear = cont_spline(s_list[0], s_list[3])
        splineshear56 = cont_spline(s56_list[0], s56_list[1])

        shearbasesum1 = (integrate(splineshear, 0, dist_s1, 1000) + integrate(splineshear, self._circumference - dist_s4, self._circumference, 1000 ))/self.skint
        shearbasesum1 -= integrate(splineshear56, -r + .001, r - 0.001, 1000)

        shearbasesum2 = (integrate(splineshear, dist_s1, dist_s1 + dist_s2, 1000) + integrate(splineshear, dist_s1 + dist_s2, dist_s1 + dist_s2 + dist_s3, 1000 ))/self.skint
        shearbasesum2 += integrate(splineshear56, -r+0.001, r - 0.001, 1000)
        #shearbasesum1 = stepcorrection * stepwidth * ((np.sum(s_list[3, idx_s[0] : idx_s[1]]) + np.sum(s_list[3, idx_s[3] : idx_s[4]]))/self.skint - np.sum(s56_list[1])/self.spart)
        #shearbasesum2 = stepcorrection *stepwidth * ((np.sum(s_list[3, idx_s[1] : idx_s[2]]) + np.sum(s_list[3, idx_s[2] : idx_s[3]]))/self.skint + np.sum(s56_list[1])/self.spart)
        twistmatrix = np.array([[2 * dist_s1 / self.skint + 2* dist_s5 / self.spart, -2* dist_s5 / self.spart ], [-2* dist_s5 / self.spart, 2 * dist_s2 / self.skint + 2 * dist_s5 / self.spart]])
        shearbasevector = np.array([[-shearbasesum1],[-shearbasesum2]])
        correctionshears = np.linalg.solve(twistmatrix, shearbasevector)
        correctionshear1 = correctionshears[0]
        correctionshear2 = correctionshears[1]


        # Add correction shears and base shear to row 4
        s_list[4, idx_s[0] : idx_s[1]] =  s_list[3, idx_s[0] : idx_s[1]] + correctionshear1
        s_list[4, idx_s[3] : idx_s[4]] =  s_list[3, idx_s[3] : idx_s[4]] + correctionshear1
        s_list[4, idx_s[1] : idx_s[2]] =  s_list[3, idx_s[1] : idx_s[2]] + correctionshear2
        s_list[4, idx_s[2] : idx_s[3]] =  s_list[3, idx_s[2] : idx_s[3]] + correctionshear2
        s56_list[4] = s56_list[1] - correctionshear1 + correctionshear2
        s_list[6] = s_list[4]
        # Arms:
        as1 = r
        as2 = np.cos( acos((self.chord - self.radius) / self.a) ) * r

        moment = as1 * stepcorrection * stepwidth *  (np.sum(s_list[3, idx_s[3]: idx_s[4]]) + np.sum(s_list[3, 0: idx_s[1]]))
        moment += np.pi*r*r/2 * correctionshear1 + r*(self.chord-r) * correctionshear2
        moment += as2 * stepcorrection * stepwidth *  (np.sum(s_list[3, idx_s[1]: idx_s[2]]) + np.sum(s_list[3, idx_s[2]: idx_s[3]+1]))
        moment = moment[0]

        # Then, the shearflow moment created about 0,0 should be equal to the moment caused by force Vy = 1 so moment = arm
        #z location wrt spar

        return moment, s_list, s56_list, correctionshears

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
            raise ValueError("The aileron does not contain this number of stiffeners")

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
        if axis == 0:
            centroid = self.centroid(2)
        else:
            centroid = self.centroid(1)
        for stiff in stiffener_list:
            d_2 = (stiff[axis] - centroid )
            i_stiff += d_2*d_2*self.stiffener_area
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
        self.xaero = np.unique(self.aerogrid[:,0])
        self._xaero = np.array([0]+list(self.xaero)[1:-1]+[self.a.span])
        self.res_locs, self.res_forces = self.get_aero_resultants(self.filename, self.aerogrid)
        self.aero_force = cont_spline(self._xaero, self.res_forces)
        self.aero_force_times_x = cont_spline(self._xaero, self.res_forces*self._xaero)
        self.res_xloc, self.totforce = self.calc_aero_res()
        self.q = partial(self._getq, self.aerogrid, self.res_forces)
        self.grid = self.make_cuts_sections(self.a, gridnum)
        self.gridsize = self.grid.size
        self.dx_list = np.asarray([self.grid[i + 1] - self.grid[i] for i in range(self.grid.size - 1)])
        self.hinge1_idx = np.argmin(np.abs(self.grid - self.a.hinge1))
        self.hinge1_val = self.grid[self.hinge1_idx]
        self.hinge2_idx = np.argmin(np.abs(self.grid - self.a.hinge2))
        self.hinge2_val = self.grid[self.hinge2_idx]
        self.hinge3_idx = np.argmin(np.abs(self.grid - self.a.hinge3))
        self.hinge3_val = self.grid[self.hinge3_idx]
        self.act1_idx = np.argmin(np.abs(self.grid - self.a.act1))
        self.act1_val = self.grid[self.act1_idx]
        self.act2_idx = np.argmin(np.abs(self.grid - self.a.act2))
        self.act2_val = self.grid[self.act2_idx]

        self.res_xloc -= self.hinge2_val

        self._geo_error = np.min(self.dx_list)*0.1

        self.hinge1z = hinge1z
        self.hinge1d_y = hinge1d_y
        self.hinge3z = hinge3z
        self.hinge3d_y = hinge3d_y
        self.P = P
        self.jammed = F_I_z
        self.hinge2_fix = hinge2_fix

        self.res_array = self.aero_sections()
        self.res_array = np.insert(self.res_array, 0, [0,0], axis=0)
        self.res_xlocs = self.res_array[:,0]
        self.res_aeroforces = self.res_array[:,1]

        # self.mac_step_vect_0 = np.vectorize(partial(mac_step, pow=0))
        # self.mac_step_vect_1 = np.vectorize(partial(mac_step, pow=1))
        self.mac_step_vect_0 = lambda difference: np.where(difference>0, 1, 0)
        self.mac_step_vect_1 = lambda difference: np.where(difference>0, difference, 0)

        if self.hinge1d_y != 0 and self.hinge3d_y != 0:
            self.bending=True
        else:
            self.bending=False

    def renew_grid(self, gridnum):
        self.grid = self.make_cuts_sections(self.a, gridnum)
        self.dx_list = [self.grid[i + 1] - self.grid[i] for i in range(self.grid.size - 1)]
        return None

    def gridnum(self, gridnum):
        self.renew_grid(gridnum)
        return None

    def grid(self, i):
        return self.grid[i]

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
        data *= 1000
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
        return np.asarray(res_locations), np.asarray(res_forces)

    def calc_aero_res(self, maxx=None):
        if maxx == None:
            maxx = self.a.span
        elif maxx < 0 or maxx > self.a.span:
            raise ValueError("Cannot calculate the resultant force for a limit outside of the domain.")
        res_force = integrate(self.aero_force, 0, maxx, int(maxx*1000))
        res_xloc = integrate(self.aero_force_times_x, 0, maxx, int(maxx*1000))/res_force
        return res_xloc, res_force

    def aero_sections(self, file="forcedata.csv"):
        res_array = np.genfromtxt(file, delimiter=",")
        return res_array

    def int_shear_y(self, x):
        return np.array([self.res_aeroforces[np.where(np.abs(self.grid - x)<self._geo_error)], - self.mac_step_vect_0(x-self.hinge1_val), -self.mac_step_vect_0(x-self.hinge2_val), -self.mac_step_vect_0(x-self.hinge3_val)])

    def int_moment_y(self, x):
        return np.array([-self.res_aeroforces[np.where(np.abs(self.grid - x)<self._geo_error)] * (x-self.res_xlocs), self.mac_step_vect_1(x-self.hinge1_val), self.mac_step_vect_1(x-self.hinge2_val), self.mac_step_vect_1(x-self.hinge3_val)])

    def int_shear_z(self, x):
        return np.array([-self.P*self.mac_step_vect_0(x-self.act2_val), self.mac_step_vect_0(x-self.hinge1_val), self.mac_step_vect_0(x-self.act1_val), self.mac_step_vect_0(x-self.hinge2_val), self.mac_step_vect_0(x-self.hinge3_val)])

    def int_moment_z(self, x):
        return np.array([self.P*self.mac_step_vect_1(x-self.act2_val), -self.mac_step_vect_1(x-self.hinge1_val), -self.mac_step_vect_1(x-self.act1_val), -self.mac_step_vect_1(x-self.hinge2_val), -self.mac_step_vect_1(x-self.hinge3_val)])

    def _getq(self, coords, forces, x):
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

    def distr_defl_func(self, q1, q2, L, I, E=73100000000):
        # deflection return positive in positive direction of q
        if q1 != q2:
            return L*L*L*L/(E*I)*(q2 /8 + (q1-q2)*3/72)
        else:
            return L*L*L*L*q2/(E*I*8)


    def distr_angle_func(self, q1, q2, L, I, E=73100000000):
        if q1 != q2:
            return L*L*L/(E*I)*(q2/6 + (q1-q2)/18)
        else:
            return L*L*L*q2/(6*E*I)

    def moment_defl(self, M, L, I, E=73100000000):
        return M*L*L*0.5/(E*I)

    def moment_angle(self, M, L, I, E=73100000000):
        return M*L/(E*I)

    def pointload_defl(self, P, L, I, E=73100000000):
        return P*L*L*L/(3*E*I)

    def pointload_angle(self, P, L, I, E=73100000000):
        return P*L*L/(2*E*I)

    def defl_point_z(self, i, force, backward):
        # force = 0 is Fz1, force = 1 is FzI, force = 2 is Fz2, force = 3 is P, force = 4 is Fz3;
        # backward = 0 is forward, backward = 1 is backward
        if force == 0:
            cutoff = self.hinge1_val
        elif force == 1:
            cutoff = self.act1_val
        elif force == 2:
            cutoff = self.hinge2_val
        elif force == 3:
            cutoff = self.act2_val
        elif force == 4:
            cutoff = self.hinge3_val
        else:
            raise ValueError("This force is not implemented for this specific function")
        if abs(self.grid[i] - cutoff) < self._geo_error:
            return self.pointload_defl(1, self.dx_list[i-backward], self.a.Iyy())
        else:
            return 0

    def angle_point_z(self, i, force, backward):
        # force = 0 is Fz1, force = 1 is FzI, force = 2 is Fz2, force = 3 is P, force = 4 is Fz3;
        # backward = 0 is forward, backward = 1 is backward
        if force == 0:
            cutoff = self.hinge1_val
        elif force == 1:
            cutoff = self.act1_val
        elif force == 2:
            cutoff = self.hinge2_val
        elif force == 3:
            cutoff = self.act2_val
        elif force == 4:
            cutoff = self.hinge3_val
        else:
            raise ValueError("This force is not implemented for this specific function")
        if abs(self.grid[i] - cutoff) < self._geo_error:
            return self.pointload_angle(1, self.dx_list[i-backward], self.a.Iyy())
        else:
            return 0

    def forward_defl_point(self, i):
        if abs(self.grid[i] - self.hinge1_val) < self._geo_error:
            return self.pointload_defl(1, self.dx_list[i], self.a.Izz())
        else:
            return 0

    def backward_defl_point(self, i):
        if abs(self.grid[i] - self.hinge3_val) < self._geo_error:
            return self.pointload_defl(1, self.dx_list[i-1], self.a.Izz())
        else:
            return 0

    def forward_angle_point(self, i):
        if abs(self.grid[i] - self.hinge1_val) < self._geo_error:
            return self.pointload_angle(1, self.dx_list[i], self.a.Izz())
        else:
            return 0

    def backward_angle_point(self, i):
        if abs(self.grid[i] - self.hinge3_val) < self._geo_error:
            return self.pointload_angle(1, self.dx_list[i-1], self.a.Izz())
        else:
            return 0

    def forward_defl_distr(self, i):
        return self.distr_defl_func(-self.q(self.grid[i + 1]), -self.q(self.grid[i]), self.dx_list[i], self.a.Izz())

    def backward_defl_distr(self, i):
        return self.distr_defl_func(-self.q(self.grid[i - 1]), -self.q(self.grid[i]), self.dx_list[i - 1], self.a.Izz())

    def forward_angle_distr(self, i):
        return self.distr_angle_func(-self.q(self.grid[i + 1]), -self.q(self.grid[i]), self.dx_list[i], self.a.Izz())

    def backward_angle_distr(self, i):
        return self.distr_angle_func(-self.q(self.grid[i - 1]), -self.q(self.grid[i]), self.dx_list[i - 1], self.a.Izz())




if __name__ == "__main__":
    a = Aileron()
    #print(a._stiff_s_pos())
    #plt.plot(a.shearcentre()[1][0], a.shearcentre()[4])
    #plt.show()
    sc = a.shearcentre(False)
    print('moment: ', sc[0])
    print('corrections: ', sc[3])
    #print("idx's: ", a.shearcentre()[2])
    plt.subplot(2, 1, 1)
    plt.plot(sc[1][0], sc[1][1], label = '1')
    plt.plot(sc[1][0], sc[1][2], label = '2')
    plt.plot(sc[1][0], sc[1][3], label = '3')
    plt.plot(sc[1][0], sc[1][4], label = '4')
    plt.subplot(2, 1, 2)
    plt.plot(sc[2][1], sc[2][0], label = '56 - 1')
    plt.plot(sc[2][4], sc[2][0], label='56 - 4')
    plt.xlabel('shearflow 56')
    plt.ylabel('y value')
    plt.legend()
    plt.show()
    #a.shearcentre()
