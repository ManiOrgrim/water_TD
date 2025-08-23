# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 07:54:19 2025

@author: manim
"""

import os
import numpy as np
import matplotlib.pyplot as plt




def readPlotScalar(filename, rocks, rockrho, poro):
    with open(filename, 'r' ) as infile:
        lines = infile.readlines()

        code = lines[2].split()[0]  # simulation code

        #read problem dimensions from line[3]
        l3ws = lines[3].split()
        Nx   = int(l3ws[0])
        Nz   = int(l3ws[1])
        Ny   = int(l3ws[2])
        Nxy  = int(l3ws[3])
        Nxz  = int(l3ws[4])
        Nxyz = int(l3ws[5])
        try:
            assert(Nx * Ny == Nxy)
            assert(Nx * Nz == Nxz)
            assert(Nx * Ny * Nz == Nxyz)
        except:
            print("occazzo")
        reachline = 4
        x = []
        y = []
        z = []
        while ("." in lines[reachline]):
            if (len(x)< Nx):
                
                x=x+[float(xvalue) for xvalue in lines[reachline].split()]
            elif(len(y)<Ny):
                y=y+[float(yvalue) for yvalue in lines[reachline].split()]
            else:
                z=z+[float(zvalue) for zvalue in lines[reachline].split()]
            reachline+=1
        
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)

        while("Time" not in lines[reachline]):  #skip boundary
            reachline +=1
        timestep = 0
        time = float(lines[reachline].split()[1])
        time_unit = (lines[reachline].split()[2])[1:-1]
        reachline +=1
        while(True):
            try:
              if ("Time" in lines[reachline]):
                 time = float(lines[reachline].split()[1])
                 time_unit = (lines[reachline].split()[2])[1:-1]
                 reachline +=1
                 timestep +=1
                 continue
              else:
                plotField(lines[reachline:reachline+int(1+Nx/12)*(Nz)+1], x, z, Nx, Nz,time,time_unit,  timestep, code, rocks, rockrho, poro)
                reachline +=int(1+Nx/12)*Nz+1 
            except (IndexError):
                break

def readPlotVector(filename, rocks):
    with open(filename, 'r' ) as infile:
        lines = infile.readlines()
 #
        code = lines[2].split()[0]  # simulation code
        #read problem dimensions from line[3]
        l3ws = lines[3].split()
        Nx   = int(l3ws[0])
        Nz   = int(l3ws[1])
        Ny   = int(l3ws[2])
        Nxy  = int(l3ws[3])
        Nxz  = int(l3ws[4])
        Nxyz = int(l3ws[5])
        try:
            assert(Nx * Ny == Nxy)
            assert(Nx * Nz == Nxz)
            assert(Nx * Ny * Nz == Nxyz)
        except:
            print("occazzo")
        reachline = 4
        x = []
        y = []
        z = []
        while ("." in lines[reachline]):
            if (len(x)< Nx):
                x=x+[float(xvalue) for xvalue in lines[reachline].split()]
            elif(len(y)<Ny):
                y=y+[float(yvalue) for yvalue in lines[reachline].split()]
            else:
                z=z+[float(zvalue) for zvalue in lines[reachline].split()]
            reachline+=1
        
 
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)

        while("Time" not in lines[reachline]):  #skip boundary
            reachline +=1
        timestep = 0
        time = float(lines[reachline].split()[1])
        time_unit = (lines[reachline].split()[2])[1:-1]
        reachline +=1
        while(True):
            try:
              if ("Time" in lines[reachline]):
                 time = float(lines[reachline].split()[1])
                 time_unit = (lines[reachline].split()[2])[1:-1]
                 reachline +=1
                 timestep +=1
                 continue
              else:
                plotVector(lines[reachline:reachline+3*int(1+Nx/12)*(Nz)+3], x, z, Nx, Nz,time,time_unit,  timestep, code, rocks)
                reachline +=3*(int(1+Nx/12)*Nz)+3 

            except (IndexError):
                break

        
def plotField (lines, x, z, Nx, Nz, time, timeunit, timestep, code, rocks, rockrho, poro):
    xfactor, xunit= convert("m", "m")
    x = x*xfactor
    zfactor, zunit= convert("m", "m")
    z = z*zfactor
    z = z[::-1]
    zorig = z[::-1]
    fieldname, fieldunit = lines[0].split(sep="(")
    fieldname = fieldname.lstrip()
    fieldunit = fieldunit.split(sep=")")[0]
    factor, fieldunit =convert(fieldunit, "MPa")
    field = []
    riga = []
    for line in lines[1:]:
       riga = riga +[float(n) for n in line.split()]
       if (len(riga)==Nx):
           field.append(riga)
           riga = []


    field = np.array(field)
    field = field*factor

    xx, zz = np.meshgrid(x,z)
    fig, ax = plt.subplots(dpi=300)
    scalefac = 1
    #◘ax.set_aspect(((max(z)-min(z))/(max(x)-min(x))))
    # fig.set_figheight(scalefac*(max(z)-min(z)))
    # fig.set_figwidth(scalefac*(max(x)-min(x)))
    if (("Fraction" in fieldname) or ("Saturation" in fieldname)):
        levels = np.linspace(0,1, 100)
    elif ("Enthalpy" in fieldname):
        levels = np.linspace(0,1000, 100)
    elif ("Pressure" in fieldname):
        levels = np.linspace(0,1.e1,100)
        
    elif ("Temper" in fieldname):
        levels = np.linspace(30, 200, 100)
    
    mappo = ax.contourf(xx, zz, field, levels=levels, extend="both")
    ax.contour(xx,zz, rocks, colors = 'black',levels=[-1.01, 0.99, 1.99, 2.99, 3.99], linewidths=2)
    ax.set_xlabel("x ["+xunit+"]")
    ax.set_ylabel("z ["+zunit+"]")
    fig.colorbar(mappo)
    ax.set_title(fieldname + " ["+fieldunit +"] at t= "+str(time) +' '+timeunit)
    ax.set_ylim(zz.max(), zz.min())
    fig.savefig(code+"_"+fieldname[:8]+'_'+"{:03d}".format(timestep)+".png")
    plt.close(fig)
    if ("Pressure" in fieldname):
        S3 = 2.180   #MPa
        Pthresh = 2*S3*(1-poro)/(3*poro*np.sqrt(poro**(-1/3)-1))
        Plith   = np.empty(field.shape)
        for k, zk in enumerate(z):
            Plith[k,:] = lithostatic_pressure(rockrho, zk)
        
        fig, ax = plt.subplots(dpi=300)
        scalefac = 1
        # ax.set_aspect(((max(z)-min(z))/(max(x)-min(x))))
        levels = np.linspace(0,1,100)
        mappo = ax.contourf(xx, zz, field/Pthresh[::-1,:], levels=np.linspace(0,1,5), extend="max")
        ax.contour(xx,zz, rocks, colors = 'black',levels=[-1.01, 0.99, 1.99, 2.99, 3.99], linewidths=2)
        ax.set_xlabel("x ["+xunit+"]")
        ax.set_ylabel("z ["+zunit+"]")
        ax.set_ylim(zz.max(), zz.min())
        fig.colorbar(mappo)
        ax.set_title("Fragmentation overpressure ratio" + " at t= "+str(time) +' '+timeunit)
        fig.savefig(code+"_"+"fragpress"+'_'+"{:03d}".format(timestep)+".png")
        plt.close(fig)
        
        
        fig, ax = plt.subplots(dpi=300)
        scalefac = 1
        # ax.set_aspect(((max(z)-min(z))/(max(x)-min(x))))
        levels = np.linspace(0,1,100)
        mappo = ax.contourf(xx, zz, field[::-1,:]/Plith[::-1,:], levels=np.linspace(0,1,5), extend="max")
        ax.contour(xx,zz, rocks, colors = 'black',levels=[-1.01, 0.99, 1.99, 2.99, 3.99], linewidths=2)
        ax.set_xlabel("x ["+xunit+"]")
        ax.set_ylabel("z ["+zunit+"]")
        ax.set_ylim(zz.max(), zz.min())
        fig.colorbar(mappo)
        ax.set_title("Lithostatic overpressure ratio" + " at t= "+str(time) +' '+timeunit)
        fig.savefig(code+"_"+"lithpress"+'_'+"{:03d}".format(timestep)+".png")
        plt.close(fig)
        
        
        
        
    
        
        
        


def plotVector (lines, x, z, Nx, Nz, time, timeunit, timestep, code, rocks):
    xfactor, xunit= convert("m", "m")
    x = x*xfactor
    zfactor, zunit= convert("m", "m")
    z = z*zfactor
    reachline= 1
    fieldname, fieldunit = lines[0].split(sep="(")
    fieldname.split("X")
    fieldname = fieldname.lstrip()
    fieldunit = fieldunit.split(sep=")")[0]
    factor, fieldunit =convert(fieldunit, fieldunit)
    fieldx = []
    riga = []
    for line in lines[reachline:reachline+int(Nx/12 +1)*Nz]:
       riga = riga +[float(n) for n in line.split()]
       if (len(riga)==Nx):
           fieldx.append(riga)
           riga = []


    fieldx = np.array(fieldx)

    reachline +=int(Nx/12+1)*Nz+1
    reachline +=int(Nx/12+1)*Nz+1
    fieldz = []
    riga = []
    for line in lines[reachline:reachline+int(Nx/12 +1)*Nz+2]:
       riga = riga +[float(n) for n in line.split()]
       if (len(riga)==Nx):
           fieldz.append(riga)
           riga = []
    fieldz = np.array(fieldz)
    xx, zz = np.meshgrid(x,z)
    fig, ax = plt.subplots(dpi=300)
    scalefac = 1
    # ax.set_aspect(((max(z)-min(z))/(max(x)-min(x))))
    # fig.set_figheight(scalefac*(max(z)-min(z)))
    # fig.set_figwidth(scalefac*(max(x)-min(x)))
    scale = 1
    # if ("Water Mass" in fieldname):
    #     scale = 1e-4
    # elif ("Water Velocity" in fieldname):
    #     scale = 0.00003
    # elif ("Steam Mass" in fieldname):
    #     scale = 0.000000001
    # elif ("Steam Velocity" in fieldname):
    #     scale = 0.00000000001
    quivero = ax.quiver(xx, zz, fieldx, fieldz, angles = "xy", scale = scale, pivot = "mid")
    ax.quiverkey(quivero, 0.8, 0.9,0.001*scale, r'{} \frac{m}{s}$', labelpos = 'E', coordinates="figure")
    ax.contour(xx,zz, rocks, colors = 'black',levels=[-1.01, 0.99, 1.99, 2.99,3.99], linewidths=1)
    ax.set_xlabel("x ["+xunit+"]")
    ax.set_ylabel("z ["+zunit+"]") 
    #fig.colorbar(mappo)
    if (timeunit== "s"):
        time = time/3600
        timeunit = "h"
    ax.set_title(fieldname + " ["+fieldunit +"] at t= "+str(time) +' '+timeunit)
    fig.savefig(code+"_"+fieldname[:8]+'_'+"{:03d}".format(timestep)+".png")
    plt.close(fig)




def convert(fieldunit, outunit):
    if ((fieldunit == "cm") and (outunit=="m")):
        factor = 0.01
    if((fieldunit=="dyne/cm^2") and (outunit=="MPa")):
        factor=1e-7
    else: 
        factor = 1
        outunit = fieldunit
        
    return factor, outunit

def ReadRock(infil, Nz):
    with open(infil, 'r') as rido:
        iline = 0
        lines = rido .readlines()
        for line in lines:
            if ("ROCK PROPERTIES") in line:
                lineadensa = lines[iline+3]
                rockdensistring = lineadensa.split()[1]
                if ',' in rockdensistring:
                    rockdensi = float(rockdensistring[:-1])
                else:
                    rocksensi = float(rockdensistring)
            if ("index_of_rock_type(i,k),i=1 to nx, k=1 to nz " in line):
                break
            else:
                iline +=1
        if ("TOP" in lines[iline+1]):
            iline+=1
        else:
            pass
        grida = []
        for line in lines[iline+1:]:
            if ("#" in line):
                break
            words = line.split()
            linea = []
            for word in words:
                if ('*' in word):
                    n, r = word.split('*')
                    for i in range(int(n)):
                        linea.append(float(r))
            grida.append(linea)
        return np.flip(np.array(grida), axis=0), rockdensi
                        
                

def read_porosity(filename):
    porosity_data = []
    inside_block = False
    
    with open(filename, "r") as f:
        for line in f:
            if "Row" in line:
                inside_block = True
                continue
            if inside_block:
                if "Row" in line:
                    continue
                if line.strip() == "":
                    # stop after the first block ends
                    break
                # collect numbers
                numbers = [float(x) for x in line.split()]
                porosity_data.extend(numbers)
    
    # Infer grid dimensions
    n_rows = int(len(porosity_data) / 16)   # since there are 16 columns (from coords seen earlier)
    porosity_array = np.array(porosity_data).reshape(n_rows, 16)
    return porosity_array


def lithostatic_pressure(rho_rock,z):
    return rho_rock*1000*9.81*z*1e-6
                
        
        
import imageio
import os
def gifferino ():
    dirac = os.listdir()
    dirac = [d for d in dirac if d[-4:]==".png"]
    images = [[dirac[0]]]
    for fname in dirac:
        if (fname[-4:] == ".png"):
            if (fname[11:19] in images[-1][-1]):           
                images[-1].append(fname)
            else:
                images.append([])
                images[-1].append(fname)
    for serie in images:
        with imageio.get_writer("{}.gif".format(serie[0][14:23]), duration=0.6, mode="I") as writer:
            for filename in serie:
                imago = imageio.imread(filename)
                writer.append_data(imago, )
        


simulaz = "JC092_00"
os.chdir("C:/Users/manim/Fare_la_Scienza/Dottorato/IAPWS/"+simulaz)
# os.chdir("t014_fold/"+simulaz)

rocks, densi = ReadRock(simulaz+".in", 5)
print("Read rock distribution")
poro=read_porosity("Out_porosity."+simulaz+".out")
readPlotScalar("Plot_scalar."+simulaz+".out", rocks, densi,  poro)
print("Plot of scalar quantities done")
# readPlotVector("Plot_vector."+simulaz+".out", rocks)
print("Plot of vector quantities done")
images = gifferino() 
print("Animated GIFs done")


        
