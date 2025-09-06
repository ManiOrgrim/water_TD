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
import re



dirs = os.listdir()
runs = [ runna for runna in dirs if ("OE" in runna and "_.npy" in runna)]

codes = [int(run[2:5]) for run in runs]

data = np.array([np.load(run) for run in runs])[:,0,1:]
i=0
# for d in data:
#     d[0]=i
#     i +=1
data = np.column_stack([codes, data])
low_lim_days =1999
data = data[data[:,-2]>= 3600*24*low_lim_days]


overlith = data[:,4].max() 
    
class Hydroplot:
    def __init__(self, data, x_index, xlabel, scale=1, logx = False, inspect =None):
        self.inspect=inspect
        self.xlabel = xlabel
        self.logx = logx
        
        self.LithData = data[:, -4]
        self.FragData = data[:, -3] 
        self.TimeData = data[:, -2]
        self.EnerData = -data[:, -1]
        self.xdata = data[:,x_index]*scale
        midf = data[data[:,-3]>0]
        
        self.LithMidf = midf[:, -4]
        self.FragMidf = midf[:, -3]
        self.TimeMidf = midf[:, -2]
        self.EnerMidf = -midf[:, -1]
        self.xmidf = midf[:,x_index]*scale
        filt = data[data[:,-3]>0.9]
        self.LithFilt = filt[:, -4]
        self.FragFilt = filt[:, -3]
        self.TimeFilt = filt[:, -2]
        self.EnerFilt = -filt[:, -1]
        self.xfilt = filt[:,x_index]*scale
        if self.inspect:
            self.LithInsp = data[self.inspect, -4]
            self.FragInsp = data[self.inspect, -3]
            self.TimeInsp = data[self.inspect, -2]
            self.EnerInsp = -data[self.inspect, -1]
            self.xinsp = data[self.inspect,x_index]*scale
    
    def plot_Frag (self, s=4):
        self.fag, self.ax = plt.subplots(dpi = 300)
        self.ax.scatter(self.xdata, self.FragData, s=s, color='limegreen')
        self.ax.scatter(self.xmidf, self.FragMidf, s=s, color='cyan')
        self.ax.scatter(self.xfilt, self.FragFilt, s=s, color='magenta')
        if self.inspect:
            self.ax.scatter(self.xinsp, self.FragInsp, s=s, color='red')
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(r"Fragmenation overpressure cells fraction")
        if (self.logx):
            self.ax.set_xscale("log")
    
    def plot_Lith (self, s=4):
        self.fag, self.ax = plt.subplots(dpi = 300)
        self.ax.scatter(self.xdata, self.LithData, s=s, color='limegreen')
        self.ax.scatter(self.xmidf, self.LithMidf, s=s, color='cyan')
        self.ax.scatter(self.xfilt, self.LithFilt, s=s, color='magenta')
        if self.inspect:
            self.ax.scatter(self.xinsp, self.LithInsp, s=s, color='red')
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(r"Lithostatic overpressure cells fraction")
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
        if self.inspect:
            self.bx1.scatter(self.xinsp, self.TimeInsp/fact, s=6, color='red')
            self.bx2.scatter(self.xinsp, self.TimeInsp/fact, s=6, color='red')
        self.bx1.spines["bottom"].set_visible(False)
        self.bx2.spines["top"].set_visible(False)
        self.bx1.xaxis.tick_top()
        self.bx1.tick_params(labeltop=False)
        self.bx2.xaxis.tick_bottom()
        self.bx2.set_ylim(0,10)
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
        if self.inspect:
            self.bx1.scatter(self.xinsp, self.EnerInsp/fact, s=6, color='red')
            self.bx2.scatter(self.xinsp, self.EnerInsp/fact, s=6, color='red')
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
        
        

# inspected = 1
        
# HeatFlux_HP = Hydroplot(data, 1, r"Heat flux [W/m$^2$]", scale=1e-3, inspect=inspected)
# HeatFlux_HP.plot_Frag()
# if overlith:
#     HeatFlux_HP.plot_Lith()
# HeatFlux_HP.plot_Time()
# HeatFlux_HP.plot_Ener()   


# Tbot_HP = Hydroplot(data, 2, r"Tbot [°C]", inspect=inspected)
# Tbot_HP.plot_Frag()
# if overlith:
#     Tbot_HP.plot_Lith()
# Tbot_HP.plot_Time()
# Tbot_HP.plot_Ener()

# TempInf_HP = Hydroplot (data, 3, r"Flux  temperature [°C]", inspect=inspected)
# TempInf_HP.plot_Frag()
# if overlith:
#     TempInf_HP.plot_Lith()
# TempInf_HP.plot_Time()
# TempInf_HP.plot_Ener()


# MassInf_HP = Hydroplot (data, 4, r"Mass influx [kg/s]", logx=True, inspect=inspected)
# MassInf_HP.plot_Frag()
# if overlith:
#     MassInf_HP.plot_Lith()
# MassInf_HP.plot_Time()
# MassInf_HP.plot_Ener()
