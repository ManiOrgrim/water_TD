import os
import subprocess
from infilewriter import modify_input_sources
import pexpect
from BOOMER import boom_eval
import numpy as np

def singlerun(name, number, flux, Tbot, HPperm, HPporo, LPperm, LPporo, inmass, intemp,  dZ):
    runname = name+"{:02d}".format(number)
    if (runname not in os.listdir()):
             os.mkdir(runname)
    os.chdir(runname)
    # with open(runname+'.in', 'w') as vraito:
    #     texto = modify_input_sources(flux, Tbot, HPperm, HPporo, LPperm, LPporo, inmass, intemp,  dZ)
    #     vraito.write(texto)
    # os.system('echo "{}" > comandi'.format(runname+'.in'))
    # os.system('echo "{}" >> comandi'.format(runname+'.out'))
    # os.system('cat "comandi" | ../ht.exe')
    print("Calculating explosivity..., run number ", name )
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
 


runs = []

#0.00000000e+00, 3.40209510e+03, 1.00092495e+03, 4.97315521e-09,
#       1.03590736e-01, 1.29044028e-17, 6.40192114e-02, 7.10886232e+00,
#       1.00000000e+00, 1.07783000e+01,            nan



simnames = ["G110C","G125C", "G150C", "G175C", "G200C", "G225C", "G250C", "G300C", "G500C", "G650C", "G800C", "G999C"]
Tbots    = [110, 125, 150, 175, 200, 225, 250, 300, 500,650, 800, 999]

# simnames = ["H150E"]
# Tbots    = [150]

def from_millidarcy2sqr_cm(millidarcy):
    return millidarcy*9.863e-12

HPperm = from_millidarcy2sqr_cm(346)
HPporo = 0.372
LPperm = 1.0e-28
LPporo = 1.0e-02#0.01+np.random.rand()*0.09
Dz =  50

# for jcond in range (len(simnames)):
#     flux = 0 #`tipycal ranges around 20000
#     runs = []
#     simname = simnames[jcond]
#     Tbot = Tbots[jcond]
#     for irun in range(41):
#         flux = irun*10000
#         OCmean, endtime, explo = singlerun(simname, irun, flux,Tbot,HPperm,HPporo,LPperm,LPporo,Dz)
#         runs.append(np.array([irun, flux,Tbot,HPperm,HPporo,LPperm,LPporo,Dz, OCmean, endtime, explo]))
#         os.chdir("..")
  
#     print(runs)
#     runs = np.array(runs)
#     np.save(simname+".npy", runs)


simnames = ["F010_100","F010_200","F010_300","F050_100","F050_200","F050_300","F100_100","F100_200","F100_300"]    
inmasses = [10,10,10,50,50,50,100,100,100]
intempes = [100,200,300,100,200,300,100,200,300]


# for jcond in range (len(simnames)):
#     flux = 0 #`tipycal ranges around 20000
#     runs = []
#     simname = simnames[jcond]
#     Tbot = 150
#     for irun in range(41):
#         flux = irun *10000
#         inmass = inmasses[jcond]
#         intemp = intempes[jcond]
#         OCmean, endtime, explo = singlerun(simname, irun, flux,Tbot,HPperm,HPporo,LPperm,LPporo,inmass, intemp, Dz)
#         runs.append(np.array([irun, flux,Tbot,HPperm,HPporo,LPperm,LPporo,Dz, OCmean, endtime, explo]))
#         os.chdir("..")
  
#     print(runs)
#     runs = np.array(runs)
#     np.save(simname+".npy", runs),    mass flux 10**[-3, 3]
seed = 2    #HC: seed 0, multisource,   mass flux 10**[-3, 3]
            #IC: seed 1, monosource     mass flux 10**[-3, 3]
            #JC: seed 2, Ttop 50°C      mass flux 10**[-3, 3]
            #KC: seed 3, Ttop 50°C      mass flux 10**[-5, 1]
np.random.seed(seed)
simchevogliovedere = [73]
for jcond in range (900):
    flux = 10000*np.random.rand()*4
    Tbot = np.random.rand()*250+50
    intemp = np.random.rand()*250+50
    inmass = 10**(np.random.rand()*6-5)
    runs = []
    simname = "JC{:03d}_".format(jcond)
    for irun in range(1):
        OCmean, endtime, explo = singlerun(simname, irun, flux,Tbot,HPperm,HPporo,LPperm,LPporo,inmass, intemp, Dz)
        runs.append(np.array([irun, flux,Tbot,intemp, inmass,HPperm,HPporo,LPperm,LPporo, Dz, OCmean, endtime, explo]))
        os.chdir("..")
  
    print(runs)
    runs = np.array(runs)
    np.save(simname+".npy", runs)
    


#singlerun("mozzd", irun, flux,Tbots,1e-10,0.1,1e-15,0.02,2)