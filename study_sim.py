# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 09:33:36 2025

@author: manim
"""

import numpy as np
import matplotlib.pyplot as plt

data = np.load("t011.npy")
data2 = np.load("tdie.npy")
data3 = np.load("t012.npy")
data4 = np.load("u012.npy")
data4m = np.load("t014.npy")


# median = np.median(np.log(data[:,-1]))
# # beans= [-1e17, -1e16, -1e15, -1e14, -1e13, -1e12, -1e10, -1e11, -1e10, 1e10, 1e11, 1e12,1e13, 1e14, 1e15, 1e16, 1e17, 1e18]
# fig, ax = plt.subplots()

# tronchi = np.log10(data[:,-1])
# tronchiN = np.log(-data[:,-1])
# ax.hist((-tronchiN))


# fEg, Ex = plt.subplots()
# Ex.hist((tronchi))


fag, ax1 = plt.subplots(dpi=300)
ax1.set_title("Explosive energy")
ax1.plot([0,1000],[1e6, 1e6], label ="Nisyros", color="green", linestyle="--")
ax1.set_ylabel("Releasable energy [MJ]")
# ax1.set_yscale("log")
ax2 = ax1.twinx()
ax2.set_ylabel("Tons of TNT")
ax2.plot([0,1000], [1800,1800], color="orange", linestyle="--", label = "Beirut 2020 ")

# ax1.plot(data[:,2], data[:,-1]*1e-6, color=[0.8,0.6, 0])
# ax2.plot(data[:,2], data[:,-1]/(4.184*1e9), color=[0.8,0.6, 0])
# ax1.plot(data2[:,2], data2[:,-1]*1e-6, color=[0.6,0.6, 0.3])
# ax2.plot(data2[:,2], data2[:,-1]/(4.184*1e9), color=[0.6,0.6, 0.3])
# ax1.plot(data3[:,2], data3[:,-1]*1e-6, color=[1,1, 0.6])
# ax2.plot(data3[:,2], data3[:,-1]/(4.184*1e9), color=[1,1, 0.6])
ax1.plot(data4m[:,2], data4m[:,-1]*1e-6, color=[0.1,0.1, 0.7])
ax2.plot(data4m[:,2], data4m[:,-1]/(4.184*1e9), color=[0.1,0.1, 0.7])


fag.legend()

