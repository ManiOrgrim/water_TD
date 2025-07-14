# -*- coding: utf-8 -*-
"""
Created on Sat Jul  5 12:46:38 2025

@author: manim
"""
import numpy as np
from sklearn.linear_model import Lasso
import iapws95 as eos
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error as MAE
from sklearn.metrics import r2_score as R2
import matplotlib.pyplot as plt

rhos_train = np.array([0.001,0.002, 0.003, 0.004, 0.005,0.006, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
                       0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.5, 1, 10, 30, 40, 50, 60, 70, 80,
                       100, 150, 160, 9 , 180, 300,321,322, 323, 232, 323, 324, 325, 600,990, 991, 992, 993, 994, 995, 997, 998,998.1, 998.4, 998.5, 998.6, 99,9,999, 999.1, 999.2, 999.3, 999.4, 999.5, 999.6,999.7,999.8,999.9, 1000,1001, 1002, 1003 ])
Ts_train   = np.array([274, 276, 279, 280, 282, 284, 290, 295, 300, 350,  400, 450, 500, 550, 600, 630, 640, 641, 642, 643, 644,645,646,647, 650, 700, 750, 900, 990, 991, 992, 993,  995, 998, 1000, 1005 ])

#rhos_train = np.concatenate((np.random.rand((1000))*100+900, np.random.rand((1000))))#, np.random.rand((1000))*900))




RHORHO, TT = np.meshgrid(rhos_train, Ts_train)
PP_t      = eos.EOS_pressure(RHORHO, TT)
rhorhoss_train = RHORHO.flatten()
TTs_train      = TT.flatten()
PPs_train      = PP_t.flatten()
ML_data=pd.DataFrame(data=zip(rhorhoss_train,TTs_train, PPs_train), columns=["rho", "T", "P"])
ML_data = ML_data[ML_data["P"]>0]
ML_data = ML_data[ML_data["P"]<1e6]
#ML_data.plot.scatter(x="T", y="P")
#ML_data.plot.scatter(x="rho", y="P")

#take from the dataframe what we want as x
x = ML_data.loc[:, ["rho", "T"]]
#take from the dataframe what we want as y
y = ML_data.loc[:, 'P']

degree = 2


p = preprocessing.PolynomialFeatures(degree).fit(x) 
print(p.get_feature_names(x.columns))
features = p.get_feature_names(x.columns)

x = preprocessing.PolynomialFeatures(degree).fit_transform(x)

x =x[:,:]
#divide x and y in train and test batches, 
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2)

#model_1 = LinearRegression().fit(x_train, y_train)
model_1 = Lasso(alpha = 0.3).fit(x_train, y_train)

y_model_train = model_1.predict(x_train)
y_model_test  = model_1.predict(x_test)

print(MAE(y_train, y_model_train))
print(MAE(y_test, y_model_test))
print(R2(y_train, y_model_train))
print(R2(y_test, y_model_test))

fig, ax = plt.subplots()
ax.plot(y_train, y_model_train, 'ro')
ax.plot(y_test, y_model_test, 'bo')
ax.plot([0,900000],[0,900000])

feg, ex = plt.subplots()
ex.bar(features, abs(model_1.coef_)) 
ex.set_yscale("log")