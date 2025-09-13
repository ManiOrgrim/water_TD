# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 09:48:17 2025

@author: manim
"""

import numpy as np
import matplotlib.pyplot as plt
from plot_ensemble import Hydroplot as HP


code ="ME"
data_frag = np.load(code+"_frag_pred.npy")
data_lith = np.load(code+"_lith_pred.npy")
#data_frag[:,0:4]
vuoto = np.zeros(data_frag[:,0].shape)
newdata = np.column_stack([data_frag[:,0:4], data_lith[:,4], data_frag[:,4], vuoto, vuoto])

# hp_heat = HP(newdata, 0, r"Heat flux [W/m$^2$]", scale=1e-3)
# hp_heat.plot_Frag(s=1)
# hp_heat.plot_Lith(s=1)
# hp_Tbot=HP(newdata, 1, r"Tbot [°C]")
# hp_Tbot.plot_Frag(s=1)
# hp_Tbot.plot_Lith(s=1)
# hp_TempInf = HP(newdata, 2, r"Flux  temperature [°C]")
# hp_TempInf.plot_Frag(s=1)
# hp_TempInf.plot_Lith(s=1)
# hp_MassInf= HP(newdata, 3, r"Mass influx [kg/s]", logx=True)
# hp_MassInf.plot_Frag(s=1)
# hp_MassInf.plot_Lith(s=1)


contdata = data_frag[::663, :]
def plot_contourf(contdatas, title, xcol, xsize, xlabel,  ycol, ysize, ylabel, xfactor=1, yfactor=1, logy=False):
    fig, ax = plt.subplots(dpi=300)
    levels = np.arange(0,1.05,0.05)
    for contdata in contdatas:
        Z = np.reshape(contdata[:,-1], (ysize,xsize))
        X = np.unique(contdata[:,xcol])/xfactor  
        Y = np.unique(contdata[:,ycol])/yfactor
        turfo = ax.contourf(Y,X,Z[::-1,:], levels=levels, alpha = 1)
        # fig.colorbar(turfo)
    ax.set_ylabel(xlabel)
    ax.set_xlabel(ylabel)
    ax.set_title(title)
    if (logy):
      ax.set_yscale("log")  
    fig.colorbar(turfo)

contdatas = [data_frag[i::663, :] for i in range(1)]
plot_contourf(contdatas, "Fragmentation overpressure", 0,51,"Heat flux [W/m$^2$]",1,51,r"T bottom [°C]",xfactor=1000)
contdata = [data_lith[660::663, :] for i in range(1)]
plot_contourf(contdatas, "Lithostatic overpressure",0,51,"Heat flux [W/m$^2$]",1,51,r"T bottom [°C]" , xfactor=1000)


contdatas = [data_frag[:663, :]]
# contdatas = [c[::13, : for c in contdatas]
plot_contourf(contdatas, "Fragmentation overpressure", 2,13, "Influx temperature [°C]",3,51,r"Influx rate [kg/s]", logy=True)
contdatas = [data_lith[:663, :]]
# contdatas = [c[::13, :] for c in contdatas]
plot_contourf(contdatas, "Lithostatic overpressure",2,13,"Influx temperature [°C]",3,51,r"Influx rate [kg/s]", logy=True )



