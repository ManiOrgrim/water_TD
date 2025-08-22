# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 13:36:03 2025

@author: manim
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib
import os


dirs = os.listdir()
runs = [ runna for runna in dirs if ("JC" in runna and "_.npy" in runna)]

data = np.array([np.load(run) for run in runs])[:,0,:]
filt = data[data[:,10]>0]

class Hydroplot:
    def __init__(self, data, x_index, xlabel, scale=1, logx = False):
        self.xlabel = xlabel
        self.logx = logx
        
        self.OverData = data[:, -3]
        self.TimeData = data[:, -2]
        self.EnerData = data[:, -1]
        self.xdata = data[:,x_index]*scale
        midf = data[data[:,10]>0]
        self.OverMidf = midf[:, -3]
        self.TimeMidf = midf[:, -2]
        self.EnerMidf = midf[:, -1]
        self.xmidf = midf[:,x_index]*scale
        filt = data[data[:,10]>0.9]
        self.OverFilt = filt[:, -3]
        self.TimeFilt = filt[:, -2]
        self.EnerFilt = filt[:, -1]
        self.xfilt = filt[:,x_index]*scale
    
    def plot_Over (self):
        self.fag, self.ax = plt.subplots(dpi = 300)
        self.ax.scatter(self.xdata, self.OverData, s=4, color='limegreen')
        self.ax.scatter(self.xmidf, self.OverMidf, s=4, color='cyan')
        self.ax.scatter(self.xfilt, self.OverFilt, s=4, color='magenta')
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(r"Overpressure cells")
        if (self.logx):
            self.ax.set_xscale("log")
    
    def plot_Time (self, timeunit='days'):
        if (timeunit == "days"):
            fact = (3600*24)
        else:
            print("Time unit not implemented")
            assert( 1==2)
        self.fbg, (self.bx1, self.bx2) = plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios':[1,3]}, dpi = 300)
        self.fbg.subplots_adjust(hspace=0.1)
        self.bx1.scatter(self.xdata, self.TimeData/fact, s=6, color='limegreen')
        self.bx2.scatter(self.xdata, self.TimeData/fact, s=6, color='limegreen')
        self.bx1.scatter(self.xmidf, self.TimeMidf/fact, s=6, color='cyan')
        self.bx2.scatter(self.xmidf, self.TimeMidf/fact, s=6, color='cyan')
        self.bx1.scatter(self.xfilt, self.TimeFilt/fact, s=6, color='magenta')
        self.bx2.scatter(self.xfilt, self.TimeFilt/fact, s=6, color='magenta')
        self.bx1.spines["bottom"].set_visible(False)
        self.bx2.spines["top"].set_visible(False)
        self.bx1.xaxis.tick_top()
        self.bx1.tick_params(labeltop=False)
        self.bx2.xaxis.tick_bottom()
        self.bx2.set_ylim(0,100)
        self.bx1.set_ylim(1999,2001)
        self.bx2.set_xlabel(self.xlabel)
        self.bx2.set_ylabel("Simulation time [{}]".format(timeunit) )
        d = 0.5
        kwargs = dict(marker=[(-1, -d), (1,d)], markersize=12, linestyle="none", mec='k', color='k', mew=1, clip_on=False)
        self.bx1.plot([0,1], [0,0], transform=self.bx1.transAxes, **kwargs)
        self.bx2.plot([0,1], [1,1], transform=self.bx2.transAxes, **kwargs)
        if (self.logx):
            self.bx1.set_xscale("log")
            self.bx2.set_xscale("log")

    def plot_Ener (self, enerunit = "MJ"):
        if (enerunit == "MJ"):
            fact = 1e6
        else:
            print("Energy unit not implemented")
            assert( 1==2)
        self.fbg, (self.bx1, self.bx2) = plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios':[0,5]}, dpi = 300)
        self.fbg.subplots_adjust(hspace=0.1)
        self.bx1.scatter(self.xdata, self.EnerData/fact, s=6, color='limegreen')
        self.bx2.scatter(self.xdata, self.EnerData/fact, s=6, color='limegreen')
        self.bx1.scatter(self.xmidf, self.EnerMidf/fact, s=6, color='cyan')
        self.bx2.scatter(self.xmidf, self.EnerMidf/fact, s=6, color='cyan')
        self.bx1.scatter(self.xfilt, self.EnerFilt/fact, s=6, color='magenta')
        self.bx2.scatter(self.xfilt, self.EnerFilt/fact, s=6, color='magenta')
        self.bx1.spines["bottom"].set_visible(False)
        self.bx2.spines["top"].set_visible(False)
        self.bx1.xaxis.tick_top()
        self.bx1.tick_params(labeltop=False)
        self.bx2.xaxis.tick_bottom()
        self.bx1.set_ylim(-2860,-2820)
        # self.bx2.set_ylim(-5000,-4200)
        self.bx2.set_xlabel(self.xlabel)
        self.bx1.set_ylabel("Releasable energy [{}]".format(enerunit) )
        d = 0.5
        kwargs = dict(marker=[(-1, -d), (1,d)], markersize=12, linestyle="none", mec='k', color='k', mew=1, clip_on=False)
        self.bx1.plot([0,1], [0,0], transform=self.bx1.transAxes, **kwargs)
        self.bx2.plot([0,1], [1,1], transform=self.bx2.transAxes, **kwargs)
        if (self.logx):
            self.bx1.set_xscale("log")
            self.bx2.set_xscale("log")
        
        
        
HeatFlux_HP = Hydroplot(data, 1, r"Heat flux [W/m$^2$]", scale=1e-3)
HeatFlux_HP.plot_Over()
HeatFlux_HP.plot_Time()
HeatFlux_HP.plot_Ener()   

Tbot_HP = Hydroplot(data, 2, r"Tbot [°C]")
Tbot_HP.plot_Over()
Tbot_HP.plot_Time()
Tbot_HP.plot_Ener()

TempInf_HP = Hydroplot (data, 3, r"Flux  temperature [°C]")
TempInf_HP.plot_Over()
TempInf_HP.plot_Time()
TempInf_HP.plot_Ener()


MassInf_HP = Hydroplot (data, 4, r"Mass influx [kg/s]", logx=True)
MassInf_HP.plot_Over()
MassInf_HP.plot_Time()
MassInf_HP.plot_Ener()
