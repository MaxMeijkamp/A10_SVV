import numpy as np
from matplotlib import pyplot as plt
import equilibrium




def comp_data(x, data_num, data_v, plot_raw=True, plot_diff_abs=True, plot_diff_rel=True, label="Verification", vertlabel="u", unit="m", fullscreen=True, extraplot=lambda: True) -> None:
    subplots=0
    if label == "Validation":
        smalllabel = "val"
    else:
        smalllabel = "ver"
    print(plot_raw)
    if plot_raw:
        rawlabel = ("Numerical model", label + " model")
        subplots += 1
        rawunit = "  ["+unit+"]"
    if plot_diff_abs:
        abslabel = ("Abs. difference numerical - "+label+" model")
        subplots += 1
        absunit = "  ["+unit+"]"
    if plot_diff_rel:
        subplots += 1
        rellabel = ("Rel. difference numerical - "+label+" model")
        relunit = "  [%]"
    if subplots == 0:
        raise ValueError("At least one plot is required for data visualization.")
    elif subplots == 1:
        if plot_raw:
            plt.plot(x, data_num)
            plt.plot(x, data_v)
            extraplot()
            plt.legend(rawlabel, loc="upper right")
            ylabel = vertlabel+rawunit
        if plot_diff_abs:
            plt.plot(x, data_num-data_v)
            plt.legend(abslabel, loc="upper right")
            plt.ylabel(vertlabel+"_num-"+vertlabel+ "_"+smalllabel+""+absunit)
        if plot_diff_rel:
            plt.plot(x, (data_num-data_v)/data_v)
            plt.legend(rellabel, loc="upper right")
            ylabel = "("+vertlabel+"_num-"+vertlabel+ "_"+smalllabel+")/("+vertlabel+"_"+smalllabel+")"+relunit
        plt.xlabel("x  [m]")
        plt.ylabel(ylabel)
        plt.show()
    elif subplots == 2:
        plot = 1
        if plot_raw:
            plt.subplot(1,2,plot)
            plt.plot(x, data_num)
            plt.plot(x, data_v)
            extraplot()
            plt.legend(rawlabel, loc="upper right")
            plt.xlabel("x  [m]")
            plt.ylabel(vertlabel+rawunit)
            plot += 1
        if plot_diff_abs:
            plt.subplot(1,2,plot)
            plt.plot(x, data_num-data_v)
            plt.legend(abslabel, loc="upper right")
            plt.xlabel("x  [m]")
            plt.ylabel(vertlabel+"_num-"+vertlabel+ "_"+smalllabel+""+absunit)
            plot += 1
        if plot_diff_rel:
            plt.subplot(1,2,plot)
            plt.plot(x, (data_num-data_v)/data_v)
            plt.legend(rellabel, loc="upper right")
            plt.xlabel("x  [m]")
            ylabel = "("+vertlabel+"_num-"+vertlabel+ "_"+smalllabel+")/("+vertlabel+"_"+smalllabel+")"+relunit
    elif subplots == 3:
        if plot_raw:
            plt.subplot(1,3,1)
            plt.plot(x, data_num)
            plt.plot(x, data_v)
            extraplot()
            plt.legend(rawlabel, loc="upper right")
            plt.xlabel("x  [m]")
            plt.ylabel(vertlabel+rawunit)
        if plot_diff_abs:
            plt.subplot(1,3,2)
            plt.plot(x, data_num-data_v)
            plt.legend(abslabel, loc="upper right")
            plt.xlabel("x  [m]")
            plt.ylabel(vertlabel+"_num-"+vertlabel+ "_"+smalllabel+""+absunit)
        if plot_diff_rel:
            plt.subplot(1,3,3)
            plt.plot(x, (data_num-data_v)/data_v)
            plt.legend(rellabel, loc="upper right")
            plt.xlabel("x  [m]")
            plt.ylabel("("+vertlabel+"_num-"+vertlabel+ "_"+smalllabel+")/("+vertlabel+"_"+smalllabel+")"+relunit)
    if fullscreen:
        figManager = plt.get_current_fig_manager()
        figManager.full_screen_toggle()
    plt.tight_layout()
    plt.show()
