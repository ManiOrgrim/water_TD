# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 14:48:27 2025

@author: manim
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib
import os
from sklearn.linear_model import LinearRegression



dirs = os.listdir()
runs = [ runna for runna in dirs if ("NE2" in runna and "_.npy" in runna)]

codes = [int(run[2:5]) for run in runs]

dataM = np.array([np.load(run) for run in runs])[:,0,1:]
i=0
dataM = np.column_stack([codes, dataM])
low_lim_days = 999
dataM = dataM[dataM[:,-2]>= 3600*24*low_lim_days]


runs = [ runna for runna in dirs if ("LB" in runna and "_.npy" in runna)]

codes = [int(run[2:5]) for run in runs]

dataL = np.array([np.load(run) for run in runs])[:,0,1:]
i=0
dataL = np.column_stack([codes, dataL])
low_lim_days =1999
dataL = dataL[dataL[:,-2]>= 3600*24*low_lim_days]

runs = [ runna for runna in dirs if ("KB" in runna and "_.npy" in runna)]
codes = [int(run[2:5]) for run in runs]

# dataK = np.array([np.load(run) for run in runs])[:,0,1:]
# i=0
# dataK = np.column_stack([codes, dataK])
# low_lim_days =1999
# dataK = dataK[dataK[:,-2]>= 3600*24*low_lim_days]



# data = np.vstack([dataK, dataM, dataL])
data = dataM


# suppose your array is called data (shape: n x 5)
X = np.column_stack([data[:, 2],data[:, 9]])   # columns x1, x2, x3, x4
# X=X.reshape(-1,1)
y = -data[:, -1] *1e-6   # column y

model = LinearRegression()
model.fit(X, y)

coefs = model.coef_
intercept = model.intercept_
r2 = model.score(X,y)

print("Coefficients:", model.coef_)
print("Intercept:", model.intercept_)
print("R² score:", model.score(X, y))

texto = r"""
E(T) = {:.2f} {:.2f} T 
R$^2$ = {:.3f}
""".format(intercept, coefs[0], r2)

fig, ax = plt.subplots(dpi=300)
ax.scatter(X[:,0], y, s=4)
ax.plot(X[:,0], model.predict(X), color='red')
ax.set_ylabel("Releasable energy [MJ]")
ax.set_xlabel("Bottom temperature [°C]")
ax.text(200, 1400, texto, fontsize=10,
    va="bottom", ha="left"  )
