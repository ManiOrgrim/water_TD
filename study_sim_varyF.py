# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 09:33:36 2025

@author: manim
"""

import numpy as np
import matplotlib.pyplot as plt


run ="A250A"
data = np.load(run+"_.npy")

# fag, ax1 = plt.subplots(dpi=300)
# ax1.set_title("Explosive energy")
# ax1.plot([0,40000],[1e6, 1e6], label ="Nisyros", color="green", linestyle="--")
# ax1.set_ylabel("Releasable energy [MJ]")
# ax2 = ax1.twinx()
# ax2.set_ylabel("Tons of TNT")
# ax2.plot([0,40000], [1800,1800], color="orange", linestyle="--", label = "Beirut 2020 ")
# ax1.plot(data[:,1], data[:,-1]*1e-6, color=[0.1,0.1, 0.7])
# fag.legend()


# feg, ex = plt.subplots(dpi=300)
# ex.set_title("Overpressure ratio")
# ex.plot(data[:,1], data[:,-3], color=[0.,0.3, 0.8], label="500 mW/m2")
# ex2 = ex.twinx()
# ex2.plot(data[:,1], data[:,-3], color='red')
# feg.legend()



# fig, ix = plt.subplots(dpi=300)
# ix.set_title("Sim time")
# ix.plot(data[1:,1], data[1:,-2]/(24*3600), color=[0.,0.3, 0.8], label="500 mW/m2")
# ix.set_ylabel("Time [d]")
# # ox.set_yscale("log")
# fig.legend()


fog, ox = plt.subplots(dpi = 300)
ox.plot (data[1:,1], data[1:,-2]/(24*3600), color="blue", label="500 mW/m2")
ox1 = ox.twinx()
ox1.plot(data[1:,1], data[1:,-1]*1e-6, color="red")
for i in range(len(data[1:,1])-1):
    ox1.axvspan((data[i,1]+data[i+1,1])/2, (data[i+1,1]+data[i+2,1])/2, color = [0,1-data[i, -3],0], alpha = 0.3)
ox.set_title (run)

