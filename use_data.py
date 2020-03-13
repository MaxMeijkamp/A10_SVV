import numpy as np
import os
from stress_modules import *
from InputClasses import Aileron

main_path = os.path.dirname(os.path.realpath(__file__))  # Script folder
data_path = os.path.join(*[main_path, 'datafiles2']) # Data folder
val_path = os.path.join(*[main_path, 'datafiles3']) # Validation folder

def dataimport(file, loadcase, data_path):
    return np.genfromtxt(data_path+"\\"+file+"_"+loadcase+".csv", delimiter=",")

class Resultsset():
    def __init__(self, loadcase, data_path, val_path, a):
        My = dataimport("My", loadcase, data_path)
        Mz = dataimport("Mz", loadcase, data_path)
        Sy = dataimport("Sy", loadcase, data_path)
        Sz = dataimport("Sz", loadcase, data_path)
        phi = dataimport("phi", loadcase, data_path)
        phi1 = dataimport("phi1", loadcase, data_path)
        phi2 = dataimport("phi2", loadcase, data_path)
        phi3 = dataimport("phi3", loadcase, data_path)
        T = dataimport("T", loadcase, data_path)
        tau = dataimport("tau", loadcase, data_path)
        v = dataimport("v", loadcase, data_path)
        v1 = dataimport("v1", loadcase, data_path)
        v2 = dataimport("v2", loadcase, data_path)
        v3 = dataimport("v3", loadcase, data_path)
        w = dataimport("w", loadcase, data_path)
        w1 = dataimport("w1", loadcase, data_path)
        w2 = dataimport("w2", loadcase, data_path)
        w3 = dataimport("w3", loadcase, data_path)
        x = np.linspace(0, a.span, int(a.span*1000+1))

        x_v = np.linspace()

        # My_v = dataimport("My", loadcase, val_path)
        # Mz_v = dataimport("Mz", loadcase, val_path)
        # Sy_v = dataimport("Sy", loadcase, val_path)
        # Sz_v = dataimport("Sz", loadcase, val_path)
        # phi_v = dataimport("phi", loadcase, val_path)
        # phi1_v = dataimport("phi1", loadcase, val_path)
        # phi2_v = dataimport("phi2", loadcase, val_path)
        # phi3_v = dataimport("phi3", loadcase, val_path)
        # T_v = dataimport("T", loadcase, val_path)
        # tau_v = dataimport("tau", loadcase, val_path)
        # v_v = dataimport("v", loadcase, val_path)
        # v1_v = dataimport("v1", loadcase, val_path)
        # v2_v = dataimport("v2", loadcase, val_path)
        # v3_v = dataimport("v3", loadcase, val_path)
        # w_v = dataimport("w", loadcase, val_path)
        # w1_v = dataimport("w1", loadcase, val_path)
        # w2_v = dataimport("w2", loadcase, val_path)
        # w3_v = dataimport("w3", loadcase, val_path)
        My_v = My
        Mz_v = Mz
        Sy_v = Sy
        Sz_v = Sz
        phi_v = phi

a = Aileron(chord = 0.605, span = 2.661, hinge1 = 0.172, hinge2 = 1.211, hinge3 = 2.591, actdist = 0.35, height = 0.205, skint = 0.0011, spart = 0.0028, stifft = 0.0012, stiffh = 0.016, stiffw = 0.019, stiffn = 15)
a.visualinspection()
for loadcase in range(1,4):
    data = Resultsset(loadcase, data_path, data_path, a)











