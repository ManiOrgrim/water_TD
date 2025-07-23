import os
import subprocess
from infilewriter import modify_input
import pexpect
from BOOMER import boom_eval
import numpy as np

def singlerun(name, number, flux, Tbot, HPperm, HPporo, LPperm, LPporo, dZ):
    runname = name+"{:04d}".format(number)
    if (runname not in os.listdir()):
             os.mkdir(runname)
    os.chdir(runname)
    with open(runname+'.in', 'w') as vraito:
        texto = modify_input(flux, Tbot, HPperm, HPporo, LPperm, LPporo, dZ)
        vraito.write(texto)
    os.system('echo "{}" > comandi'.format(runname+'.in'))
    os.system('echo "{}" >> comandi'.format(runname+'.out'))
    os.system('cat "comandi" | ../ht.exe')
    print("Calculating explosivity..., run number ", number )
    overcome, runtime, explosivity = boom_eval(runname)
    #print(explosivity)
    return overcome, runtime, -explosivity[0]

    


fluxes = [0, 1, 10, 50, 100, ] #B-ok, T
Tbots  = [ 30, 50, 100, 600  ]
HPperms= [ 1e-9, 1e-10, 1e-11]
HPporos= [ 0.1, 0.15]
LPperms= [ 1e-17, 1e-19]
LPporos= [0.05, 0.01]
dZs    = [ 5, 10, 15, 20]


  

#fluxes = [0, 10, 100, 500, 1000] 
#Tbots  = [ 30, 10, 400, 600,900 ]
#HPperms= [ 1e-10, 1e-9, 1e-8]
#HPporos= [ 0.15, 0.2]
#LPperms= [1e-12, 1e-17, 1e-19, 1e-20]
#LPporos= [0.05, 0.01]
#dZs    = [ 1,3, 5, 7, 10, 20]
 


simname= "t013"
runs = []

#0.00000000e+00, 3.40209510e+03, 1.00092495e+03, 4.97315521e-09,
#       1.03590736e-01, 1.29044028e-17, 6.40192114e-02, 7.10886232e+00,
#       1.00000000e+00, 1.07783000e+01,            nan



HPperm = 1e-14
HPporo = 0.15
LPperm = 1.0e-19
LPporo = 5.0e-02#0.01+np.random.rand()*0.09
Dz =  10
Tbot = 30
flux = 1000


#for irun in range(3000):
for irun in range(1000):
    # flux = 3.40209510e+03#np.random.rand()*4000
    # # Tbot = 1.00092495e+03#30+np.random.rand()+970
    # HPperm = 1e-14
    # HPporo = 0.15
    # LPperm = 1.0e-19
    # LPporo = 5.0e-02#0.01+np.random.rand()*0.09
    # Dz =  10
    Tbot +=1
    # HPporo = 0.1+np.random.rand()*0.15
    # LPperm = 10**(-12-np.random.rand()*8)
    #LPporo = 0.01+np.random.rand()*0.09
    # Dz =  1+np.random.rand()*29

    
    OCmean, endtime, explo = singlerun(simname, irun, flux,Tbot,HPperm,HPporo,LPperm,LPporo,Dz)
    runs.append(np.array([irun, flux,Tbot,HPperm,HPporo,LPperm,LPporo,Dz, OCmean, endtime, explo]))
    os.chdir("..")
  
print(runs)
runs = np.array(runs)
np.save(simname+".npy", runs)

#singlerun("mozzd", irun, flux,Tbots,1e-10,0.1,1e-15,0.02,2)