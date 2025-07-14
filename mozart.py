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
    print("Calculating explosivity...")
    overcome, runtime, explosivity = boom_eval(runname)
    print(explosivity)
    return overcome, runtime, explosivity[0]

    


fluxes = [0, 1, 10, 50, 100, ] #B-ok, T
Tbots  = [ 30, 50, 100, 600  ]
HPperms= [ 1e-9, 1e-10, 1e-11]
HPporos= [ 0.1, 0.15]
LPperms= [ 1e-17, 1e-19]
LPporos= [0.05, 0.01]
dZs    = [ 5, 10, 15, 20]


  

fluxes = [0, 10, 100, 500, 1000] 
Tbots  = [ 30, 10, 400, 600,900 ]
HPperms= [ 1e-10, 1e-9, 1e-8]
HPporos= [ 0.15, 0.2]
LPperms= [1e-12, 1e-17, 1e-19, 1e-20]
LPporos= [0.05, 0.01]
dZs    = [ 1,3, 5, 7, 10, 20]
 


simname= "tqat"
irun = 0
runs = []
for flux in fluxes:
    for Tbot in Tbots:
        for HPperm in HPperms:
            for HPporo in HPporos:
                for LPperm in LPperms:
                    for LPporo in LPporos:
                        for dZ in dZs:
                            OCmean, endtime, explo = singlerun(simname, irun, flux,Tbot,HPperm,HPporo,LPperm,LPporo,dZ)
                            runs.append(np.array([irun, flux,Tbot,HPperm,HPporo,LPperm,LPporo,dZ, OCmean, endtime, explo]))
                            print(runs)
                            os.chdir("..")
  
                            irun +=1
print(runs)
runs = np.array(runs)
np.save(simname+".npy", runs)

#singlerun("mozzd", irun, flux,Tbots,1e-10,0.1,1e-15,0.02,2)