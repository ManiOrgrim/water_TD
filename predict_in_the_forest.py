# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 07:16:58 2025

@author: manim
"""

import numpy as np
import os
from collections import Counter
from sklearn import datasets
from sklearn.model_selection import train_test_split
from RandomForest import RandomForest

dirs = os.listdir()
runs = [ runna for runna in dirs if ("KB" in runna and "_.npy" in runna)]

data = np.array([np.load(run) for run in runs])[:,0,:]
low_lim_days =1999
data = data[data[:,-2]>= 3600*24*low_lim_days]

X, y = data[:,1:5], data[:,10]
X_train, X_test, y_train, y_test = train_test_split(X,y, test_size=0.2, random_state=1234)

print("Planting the forest")
clf = RandomForest()

clf.fit(X_train, y_train)
# predictions, stds = clf.predict(X_test)
print("Fit done")
def accuracy(y_test, y_pred):
    return np.sqrt(np.sum((y_test- y_pred)**2))/np.sum(y_test)


# predictions, stds = clf.predict(X_test)

fluxes = np.linspace(0, 40000, 51) # 0- 40_000
Tbots = np.linspace(50, 300, 51)#50-300
intemps = np.linspace(50, 300, 51)#50-300
inmasses = 10**(np.linspace(-5,1,13))#10**[-5, 1]


X = []

print("Building the input array")
i = 0
for flux in fluxes:
    i+=1
    for Tbot in Tbots:
        for intemp in intemps:
            for inmass in inmasses:
                X.append(np.array([flux, Tbot, intemp, inmass]))

print("Predicting")
pred,_ = clf.predict(X)
sets = np.column_stack([X, pred])
np.save("KB_lith_pred.npy", sets)