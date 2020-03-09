import numpy as np
from numericaltools import *
from InputClasses import *
from equilibrium import *
from compare import *
import os
import sys
from matplotlib import rc


# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)


################# Toggle functionality ##################
# Numerical model data generation
renew_data_y = False
renew_data_z = False
# Verification model data generation
verification_gen = False
# Validation model data generation
validation_gen = False

# Toggle comparisons
verify = True
validate = True

main_path = os.path.dirname(os.path.realpath(__file__))  # Script folder
data_path = os.path.join(*[main_path, 'datafiles']) # Data folder

a = Aileron()
loads = AppliedLoads(int(a.span*1000), aileron=a)

if renew_data_y:
    uy, sheary, momenty = make_u_y(loads)
    uy_new, xy = calc_deflection_y(loads, uy)
    shearsolvedy = sheary[:,0] + sheary[:, 1:].dot(xy[1:])
    momentsolvedy = momenty[:,0] + momenty[:, 1:].dot(xy[1:])
    uy_save = np.copy(uy)
    uy_save[:,1] *= xy[0]
    uy_save[:,2] *= xy[1]
    uy_save[:,3] *= xy[2]
    uy_save[:,4] *= xy[3]
    save_data(xy, appendix="y", prefix=str(data_path)+"\\", shear_blank=sheary, shear_solved=shearsolvedy,
              moment_blank=momenty, moment_solved=momentsolvedy, u_blank=uy, u=uy_save)

if renew_data_z:
    uz, shearz, momentz = make_u_z(loads)
    uz_new, xz = calc_deflection_z(loads, uz)
    shearsolvedz = shearz[:,0] + shearz[:, 1:].dot(xz[1:])
    momentsolvedz = momentz[:,0] + momentz[:, 1:].dot(xz[1:])
    uz_save = np.copy(uz)
    uz_save[:,1] *= xz[0]
    uz_save[:,2] *= xz[1]
    uz_save[:,3] *= xz[2]
    uz_save[:,4] *= xz[3]
    uz_save[:,5] *= xz[4]
    save_data(xz, appendix="z", prefix=str(data_path)+"\\", shear_blank=shearz, shear_solved=shearsolvedz,
              moment_blank=momentz, moment_solved=momentsolvedz, u_blank=uz, u=uz_save)

if verification_gen: # Does not work 100% yet
    if sys.version_info.minor == 7:
        from A10_SVV_VerificationModel.main import *
    elif sys.version_info.minor == 6:
        from A10_SVV_VerificationModel_3_6.main import *
    else:
        raise ValueError("This python version's verification model has not been downloaded.")

    ####### Following code thanks to the verification model: ########
    # van Elsloo, S.J., van Campen, J.M.J.F & van der Wal, W. (2020). Verification model
    x = np.linspace(0,a.span, int(a.span*1000 + 1))
    v, w, phi = aileron.eval(x)       # Compute the three deflections
    v1, w1, phi1 = aileron.fdeval(x)     # Compute their their first order derivative
    v2, w2, phi2 = aileron.sdeval(x)     # Compute their their second order derivative
    v3, w3, phi3 = aileron.tdeval(x)     # Compute their their third order derivative
    # Compute the loading
    Sy = aileron.Sy(x)  # Compute the shear force in y
    Sz = aileron.Sz(x)  # Compute the shear force in z
    My = aileron.My(x)  # Compute the moment around the y-axis
    Mz = aileron.Mz(x)  # Compute the moment around the z-axis
    T = aileron.T(x)  # Compute the torque
    tau = aileron.tau(x)  # Compute the distributed torque

    np.savetxt(str(data_path)+"\\S_y.csv", Sy, delimiter=",")
    np.savetxt(str(data_path)+"\\S_z.csv", Sy, delimiter=",")
    np.savetxt(str(data_path)+"\\M_y.csv", Sy, delimiter=",")
    np.savetxt(str(data_path)+"\\M_z.csv", Sy, delimiter=",")
    np.savetxt(str(data_path)+"\\T.csv", Sy, delimiter=",")
    np.savetxt(str(data_path)+"\\tau.csv", Sy, delimiter=",")

def dataimport(file):
    return np.genfromtxt(data_path+"/"+file+".csv", delimiter=",")

sol_y = dataimport("x_y") # Angle at D wrt x axis, Fy1, Fy2, Fy3
sol_z = dataimport("x_z") # Angle at D wrt z axis, Fz1, Fza, Fz2, Fz3
x = np.linspace(0, a.span, int(a.span*1000+1))

Sy_num = dataimport("shearsolved_y")
Sy_ver = dataimport("Sy")
My_num = dataimport("momentsolved_y")
My_ver = dataimport("My")
Mz_num = dataimport("momentsolved_z")
Mz_ver = dataimport("Mz")
uy_num = np.sum(dataimport("usolved_y"),axis=1)
uy_ver = dataimport("v")
uz_num = np.sum(dataimport("usolved_z"),axis=1)
uz_ver = dataimport("usolved_z")

# uz_val = dataimport("VM_Jam_bent_values")
# print(uz_val)
# b = dataimport("VM_Jam_straight_coordinates")[:,0]
# func = cont_spline(x, uz_num)
# values = func(b)
# comp_data(b, uz_num, uz_val, plot_diff=False, plot_rel=False, label="Validation")

uy1_num = derivative(uy_num, x)
uy1_ver = dataimport("v1")
uy2_num = derivative(uy1_num, x)
uy2_ver = dataimport("v2")
uy3_num = derivative(uy2_num, x)
uy3_ver = dataimport("v3")

def plot_bc():
    plt.scatter(loads.hinge1_val, loads.hinge1d_y)
    plt.scatter(loads.hinge2_val, 0)
    plt.scatter(loads.hinge3_val, loads.hinge3d_y)

# comp_data(x, Sy_num, Sy_ver, vertlabel="Sy", unit="N")
# comp_data(x, My_num, My_ver, vertlabel="My", unit="Nm")
# comp_data(x, uz_num, uz_ver, vertlabel="u", unit="m")
comp_data(x, uy3_num, uy3_ver, vertlabel="u", unit="m")#, extraplot=plot_bc)
