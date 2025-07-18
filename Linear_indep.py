# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 09:29:29 2025

@author: manim
"""

import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Lasso
from sklearn.metrics import mean_absolute_error as MAE
from sklearn.metrics import r2_score as R2
import matplotlib.pyplot as plt
from sklearn.preprocessing import FunctionTransformer



data = np.load("tott.npy")

ML_data  = pd.DataFrame(data=data, columns = ["Nsim", "Heat_Flux", "Tbot", "HPperm", "HPporo", "LPperm", "LPporo", "Dz", "Overpressure_perc", "Endtime", "Explo_energy" ])
ML_data["Thermonum"]   = ((10**4) * ML_data["Heat_Flux"]*np.sqrt(ML_data["HPperm"])/(1.5e5 *(ML_data["Tbot"]-30)))
ML_data["Thermonum2"]  = ((10**4) * ML_data["Heat_Flux"]*np.sqrt(ML_data["LPperm"])/(1.5e5 *(ML_data["Tbot"]-30)))**(-1)
ML_data["Thermonum3"]  = ((10**4) * ML_data["Heat_Flux"]*ML_data["Dz"]/(1.5e5 *(ML_data["Tbot"]-30)))
# ML_data["Thermonum"]   = (ML_data["HPporo"]*(10**4) * ML_data["Heat_Flux"]*np.sqrt(ML_data["HPperm"])/(1.5e5 *(ML_data["Tbot"]-30)))**(-1)
# ML_data["Thermonum2"]  = (ML_data["HPporo"]*(10**4) * ML_data["Heat_Flux"]*np.sqrt(ML_data["LPperm"])/(1.5e5 *(ML_data["Tbot"]-30)))**(-1)

ML_data["log_HPperm"]  = np.log10(ML_data["HPperm"])
ML_data["log_LPperm"]  = np.log10(ML_data["LPperm"])
ML_data["perme_ratio"] = ML_data["HPperm"]/ML_data["LPperm"]


# take from the dataframe the inputs

x_orig = ML_data.loc[:, ["Heat_Flux", "Tbot", "log_HPperm", "HPporo", "log_LPperm", "LPporo", "Dz", "Thermonum", "Thermonum2"]]
x = np.zeros(x_orig.shape)
x=x_orig


# logs of permeabilities are actually calculated here
# x["log_HPperm"] = np.log10(x["log_HPperm"])
# x["log_LPperm"] = np.log10(x["log_LPperm"])




# adimensionalizers = np.array([0.5*4000, 30+0.5*970, np.log10(10**(-7-0.5*3)), 0.1+0.5*0.15, np.log10(10**(-12-0.5*8)), 0.01+0.5*0.09, 1+0.5*29])
# x = x[:]/adimensionalizers


# start working on the energy
y = ML_data.loc[:, "Explo_energy"]

degree = 4


p = preprocessing.PolynomialFeatures(degree).fit(x)

features = p.get_feature_names(x.columns)
scaler = preprocessing.StandardScaler().fit(x)

x = scaler.transform(x) 
x = preprocessing.PolynomialFeatures(degree).fit_transform(x)


# MODEL 1
x_train, x_test, y_train, y_test = train_test_split(x,y,test_size=0.2)

model_01 = Lasso(alpha = 0.0).fit(x_train, y_train)

y_model_train = model_01.predict(x_train)
y_model_test  = model_01.predict(x_test)

print("_______________________________________________________________")
print("model_01")
print("MAE train", MAE(y_train, y_model_train))
print("MAE test", MAE(y_test, y_model_test))
print("R2 train", R2(y_train, y_model_train))
print("R2 test", R2(y_test, y_model_test))


fagmodel_01, axmodel_01 = plt.subplots(2)
axmodel_01[0].plot(y_train, y_model_train, 'ro')
axmodel_01[0].plot(y_test, y_model_test, 'bo')
axmodel_01[0].plot([0,900000],[0,900000])
axmodel_01[1].bar(features, abs(model_01.coef_))
axmodel_01[0].set_title("model_01") 


minidx = np.argmin(abs(model_01.coef_ )) 

print("Minimal coefficient is coeff number", minidx )
print(features[minidx], "with values", model_01.coef_ [minidx])





# MODEL 2 - remove feature  

for i in range(810):
    if (i> 0):
       x  = np.delete(x, minidx, axis=1)
       features.pop(minidx)
    
    x_train, x_test, y_train, y_test = train_test_split(x,y,test_size=0.2)
    model_02 = Lasso(alpha = 0.1+0.07*i).fit(x_train, y_train)
    
    
    
    y_model_train = model_02.predict(x_train)
    y_model_test  = model_02.predict(x_test)
    print("_______________________________________________________________")
    print("model ", i)
    print("MAE train", MAE(y_train, y_model_train))
    print("MAE test", MAE(y_test, y_model_test))
    print("R2 train", R2(y_train, y_model_train))
    print("R2 test", R2(y_test, y_model_test))
    if ( R2(y_train, y_model_train)<0.75):
        fagmodel_02, axmodel_02 = plt.subplots(2)
        axmodel_02[0].plot(y_train, y_model_train, 'ro')
        axmodel_02[0].plot(y_test, y_model_test, 'bo')
        axmodel_02[0].plot([0,900000],[0,900000])
        axmodel_02[1].bar(features, abs(model_02.coef_))
        axmodel_02[0].set_title("model_{}".format(i))
        break
    

    abso = abs(model_02.coef_)
    minidx = np.argmin(abso)
    
    print("Minimal coefficient is coeff number", minidx )
    print(features[minidx], "with values {:.2e}".format( model_02.coef_ [minidx]))







