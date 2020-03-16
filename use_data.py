import numpy as np
import os
from stress_modules import *
from InputClasses import Aileron
from compare import *
from math import asin, acos

main_path = os.path.dirname(os.path.realpath(__file__))  # Script folder
data_path = os.path.join(*[main_path, 'datafiles2']) # Data folder
val_path = os.path.join(*[main_path, 'datafileval']) # Validation folder


def dataimport(file, loadcase, data_path):
    return np.genfromtxt(data_path+"\\"+file+"_"+loadcase+".csv", delimiter=",")

def changex(xold, dataold, xnew):
    temp = cont_spline(xold, dataold)
    return temp(xnew)

def calc_s(a: Aileron, y, z):
    if z>-a.radius+0.000001:
        if y>=0:
            return a.radius * asin(y/a.radius), False
        else:
            return a._circumference + a.radius * asin(y/a.radius), False
    if z<-a.radius-0.000001:
        if y>=0:
            return a.radius*np.pi*0.5 + (a.radius-y)*a.a/a.radius, False
        if y<0:
            return a._circumference*0.5-y*a.a/a.radius , False
    if abs(z+a.radius)<0.000001:
        return y+a.radius, True

def calc_yz(a, s):
    if s <= a.radius*np.pi*0.5:
        return a.radius*sin(s/a.radius), a.radius*cos(s/a.radius)
    elif s >= a._circumference-a.radius*np.pi*0.5:
        return -a.radius*sin((a._circumference-s)/a.radius), a.radius*cos((a._circumference-s)/a.radius)
    elif s > a._circumference*0.5:
        return -(s-a._circumference*0.5)*a.radius/a.a, (s-a._circumference*0.5)*(a.chord-a.radius)/a.a-a.chord+a.radius
    elif s > a.radius*np.pi*0.5:
        s_part = s-a.radius*np.pi*0.5
        return a.radius-(s_part*a.radius/a.a), -(s_part*(a.chord-a.radius)/a.a)


calc_s_vect = np.vectorize(calc_s)
calc_yz_vect = np.vectorize(calc_yz)


def calc_yz_total(a, s, sparidx=-100):
    result1 = calc_yz_vect(a, s[:sparidx])
    result2 = s[sparidx:]-a.radius, np.zeros(s[sparidx:].size)
    return np.hstack((result1, result2))[::-1]


class Resultsset():
    def __init__(self, loadcase, data_path, val_path, a):
        self.My = dataimport("My", loadcase, data_path)
        self.Mz = dataimport("Mz", loadcase, data_path)
        self.Sy = dataimport("Sy", loadcase, data_path)
        self.Sz = dataimport("Sz", loadcase, data_path)
        self.phi = dataimport("phi", loadcase, data_path)
        self.phi1 = dataimport("phi1", loadcase, data_path)
        self.phi2 = dataimport("phi2", loadcase, data_path)
        self.phi3 = dataimport("phi3", loadcase, data_path)
        self.T = dataimport("T", loadcase, data_path)
        self.tau = dataimport("tau", loadcase, data_path)
        self.v = dataimport("v", loadcase, data_path)
        self.v1 = dataimport("v1", loadcase, data_path)
        self.v2 = dataimport("v2", loadcase, data_path)
        self.v3 = dataimport("v3", loadcase, data_path)
        self.w = dataimport("w", loadcase, data_path)
        self.w1 = dataimport("w1", loadcase, data_path)
        self.w2 = dataimport("w2", loadcase, data_path)
        self.w3 = dataimport("w3", loadcase, data_path)

        self.s = np.linspace(0, a._circumference, 1000)
        self.s_alt = np.linspace(0, a.height, 100)

        self.c1 = dataimport("c1", loadcase, data_path)
        z1 = dataimport("d1", loadcase, data_path)
        y1 = dataimport("e1", loadcase, data_path)
        self.s1 = a.radius*np.linspace(0,np.pi/2,num = 100)
        self.c2 = dataimport("c2", loadcase, data_path)
        z2 = dataimport("d2", loadcase, data_path)
        y2 = dataimport("e2", loadcase, data_path)
        self.s2, _ = calc_s_vect(a, y2-a.radius, z2)
        self.c3 = dataimport("c3", loadcase, data_path)
        z3 = dataimport("d3", loadcase, data_path)
        y3 = dataimport("e3", loadcase, data_path)
        self.s3 = np.linspace(np.pi*0.5*a.radius,np.pi*0.5*a.radius+a.a,num = 100)
        self.c4 = dataimport("c4", loadcase, data_path)
        z4 = dataimport("d4", loadcase, data_path)
        y4 = dataimport("e4", loadcase, data_path)
        self.s4 = np.linspace(a._circumference*0.5 ,a._circumference*0.5+a.a,num = 100)
        self.c5 = dataimport("c5", loadcase, data_path)
        z5 = dataimport("d5", loadcase, data_path)
        y5 = dataimport("e5", loadcase, data_path)
        self.s5, _ = calc_s_vect(a, y5, z5)
        self.s5 -= a.radius
        self.c6 = dataimport("c6", loadcase, data_path)
        z6 = dataimport("d6", loadcase, data_path)
        y6 = dataimport("e6", loadcase, data_path)
        self.s6 = a.radius*np.linspace(-np.pi/2,0,num = 100)+a._circumference
        self.c_skin = np.concatenate([self.c1, self.c3[1:], self.c4[1:], self.c6[1:]])
        self.s_skin = np.concatenate([self.s1, self.s3[1:], self.s4[1:], self.s6[1:]])
        self.c_skin = self.c_skin[np.argsort(self.s_skin)]
        self.s_skin = self.s_skin[np.argsort(self.s_skin)]
        self.c_skin = changex(self.s_skin, self.c_skin, self.s)
        self.s_spar = np.concatenate([self.s5, self.s2[1:]])
        self.c_spar = np.concatenate([self.c5, self.c2[1:]])
        self.c_spar = self.c_spar[np.argsort(self.s_spar)][::-1]
        self.s_spar = self.s_spar[np.argsort(self.s_spar)]
        self.c_spar = changex(self.s_spar, self.c_spar, self.s_alt)

        zsc = -0.10856995078063854
        self.vh = self.v
        self.wh = self.w + self.phi * (zsc+a.height/2)

        self.x = np.linspace(0, a.span, int(a.span*1000+1))

        self.x_v = np.genfromtxt(val_path+"\\Xcoord")
        self.x_v = self.x_v*0.001

        phi_v = np.genfromtxt(val_path+"\\phi_"+loadcase)
        self.phi_v = changex(self.x_v, phi_v, self.x)
        v_v = np.genfromtxt(val_path+"\\v_"+loadcase)
        self.v_v = changex(self.x_v, v_v, self.x)*0.001
        w_v = np.genfromtxt(val_path+"\\w_"+loadcase)
        self.w_v = changex(self.x_v, w_v, self.x)*0.001

        VMc = np.genfromtxt(val_path+"\\VMc_"+loadcase+".csv", delimiter=",")
        self.vmy, self.vmz = VMc[:,1]*0.001, VMc[:,2]*0.001-a.radius

        self.vms, self.spar = calc_s_vect(a, self.vmy, self.vmz)
        self.vm = np.genfromtxt(val_path+"\\VMv_"+loadcase+".csv", delimiter=",")
        self.vmskin = self.vm[~self.spar]
        self.vmspar = self.vm[self.spar]
        self.vms_skin = self.vms[~self.spar]
        self.vms_spar = self.vms[self.spar]
        self.vmskin = self.vmskin[np.argsort(self.vms_skin)]
        self.vms_skin = self.vms_skin[np.argsort(self.vms_skin)]
        self.vmspar = self.vmspar[np.argsort(self.vms_spar)]
        self.vms_spar = self.vms_spar[np.argsort(self.vms_spar)]
        extra0 = (self.vmskin[0]+self.vmskin[-1])*0.5
        self.vmskin = np.concatenate((np.array([extra0]), self.vmskin, np.array([extra0])))
        self.vms_skin = np.concatenate((np.array([0]), self.vms_skin, np.array([a._circumference])))
        extra0 = self.vmspar[0]-(self.vmspar[1]-self.vmspar[0])/(self.vms_spar[1]-self.vms_spar[0])*self.vms_spar[0]
        extra1 = self.vmspar[-1]+(self.vmspar[-1]-self.vmspar[-2])/(self.vms_spar[-1]-self.vms_spar[-2])*(a.height-self.vms_spar[-1])
        self.vmspar = np.concatenate((np.array([extra0]), self.vmspar, np.array([extra1])))
        self.vms_spar = np.concatenate((np.array([0]), self.vms_spar, np.array([a.height])))
        self.vmskin = changex(self.vms_skin, self.vmskin, self.s)*1e9
        self.vmspar = changex(self.vms_spar, self.vmspar, self.s_alt)*1e9
        self.vdata_val = np.concatenate((self.vmskin, self.vmspar))
        self.vdata_ver = np.concatenate((self.c_skin, self.c_spar))
        self.sdata = np.concatenate((self.s, self.s_alt))


a = Aileron(chord = 0.605, span = 2.661, hinge1 = 0.172, hinge2 = 1.211, hinge3 = 2.591, actdist = 0.35, height = 0.205, skint = 0.0011, spart = 0.0028, stifft = 0.0012, stiffh = 0.016, stiffw = 0.019, stiffn = 15)

list1 = []
list2 = []

for loadcase in range(1,4):
    print("Loadcase =", loadcase)
    data = Resultsset(str(loadcase), data_path, val_path, a)
    # z, y = calc_yz_total(a, data.sdata, sparidx=-100)
    # plt.subplot(221)
    # plt.scatter(z, y, c=(data.vdata_ver))
    # plt.colorbar()
    # plt.axis('equal')
    # plt.title("Verification data")
    # plt.subplot(222)
    # plt.scatter(z, y, c=(data.vdata_val))
    # plt.colorbar()
    # plt.axis('equal')
    # plt.title("Validation data")
    # plt.subplot(223)
    # plt.scatter(z, y, c=(data.vdata_ver-data.vdata_val))
    # plt.colorbar()
    # plt.axis('equal')
    # plt.title("Absolute difference")
    # plt.subplot(224)
    # plt.scatter(z, y, c=(data.vdata_ver-data.vdata_val)/data.vdata_val)
    # plt.colorbar()
    # plt.axis('equal')
    # plt.title("Relative difference")
    # figManager = plt.get_current_fig_manager()
    # figManager.full_screen_toggle()
    # plt.tight_layout()
    # plt.show()
    comp_data(data.x, data.phi-np.average(data.phi), -data.phi_v+np.average(data.phi_v), label="Validation", label0="Verification", vertlabel="phi", unit="rad", fullscreen=False)
    # comp_data(data.x, data.vh, data.v_v, label="Validation", label0="Verification", vertlabel="v", unit="m", fullscreen=False)
    # comp_data(data.x, data.w, data.w_v, label="Validation", label0="Verification", vertlabel="w", unit="m", fullscreen=False)












