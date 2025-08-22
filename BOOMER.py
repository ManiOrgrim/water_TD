# -*- coding: utf-8 -*-
"""
Created on Sun Jul  6 08:23:31 2025

@author: manim
"""

import iapws95 as eos

import numpy as np

def find_last_occurrence(filepath, search_string):
    """
    Returns the line index (0-based) of the last line in the file that contains search_string.
    Returns None if the string is not found.
    """
    last_index = None
    with open(filepath, "r") as file:
        for i, line in enumerate(file):
            if search_string in line:
                last_index = i
    return last_index



def read_temperature_file(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()

    # === Extract X Coordinates ===
    x_start = lines.index(next(l for l in lines if "R-Direction Node Coordinate" in l)) + 3
    x_coords = []
    while lines[x_start].strip():
        parts = lines[x_start].split()
        x_coords.extend([float(p) for p in parts])
        x_start += 1
    x_start = lines.index(next(l for l in lines if "R-Direction Node Coordinate" in l)) + 8
    while lines[x_start].strip():
        parts = lines[x_start].split()
        x_coords.extend([float(p) for p in parts])
        x_start += 1

    # === Extract Y Coordinates ===
    try:
      y_start = lines.index(next(l for l in lines if "Y-Direction Node Coordinate" in l)) + 2
      y_coords = []
      while lines[y_start].strip():
        parts = lines[y_start].split()
        y_coords.extend([float(p) for p in parts])
        y_start += 1
    except:
        y_coords =[1]
    # === Extract Z Coordinates ===
    z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 3
    z_coords = []
    while lines[z_start].strip():
        parts = lines[z_start].split()
        z_coords.extend([float(p) for p in parts])
        z_start += 1
    z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 8
    while lines[z_start].strip():
        parts = lines[z_start].split()
        z_coords.extend([float(p) for p in parts])
        z_start += 1
    # z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 13
    # while lines[z_start].strip():
    #     parts = lines[z_start].split()
    #     z_coords.extend([float(p) for p in parts])
    #     z_start += 1
    # z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 18
    # while lines[z_start].strip():
    #     parts = lines[z_start].split()
    #     z_coords.extend([float(p) for p in parts])
    #     z_start += 1

    # === Extract Temperature Field ===
    #temp_start = lines.index(next(l for l in lines[::-1] if "--- Temperature Values" in l)) + 2
    temp_start =find_last_occurrence(filepath, "Temperature")+3
    temperature = []
    for line in lines[temp_start:]:
        parts = line.strip().split()
        if len(parts) < 2:
            break
        row = [float(p) for p in parts[1:]]  # skip Z-index!exitr
        temperature.append(row)
    for i in range(len(temperature)):
        line = lines[temp_start+4+len(z_coords)+i]
        parts = line.strip().split()
        if len(parts) < 2:
            break
        temperature[i]=np.array(temperature[i]+[float(p) for p in parts[1:]])  # skip Z-index))
        # temperature.append(row)

    # Convert everything to numpy arrays
    x = np.array(x_coords)
    y = np.array(y_coords)
    z = np.array(z_coords)
    T = np.array(temperature)  # shape should be (len(z), len(x))

    return x, y, z, T 

def read_porosity(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    poro_start =find_last_occurrence(filepath, "Media Porosity")+2
    poro = []
    poroline = []
    for line in lines[poro_start:]:
        if ("Row") in line:
            if (len(poroline)>0):
               poro.append(poroline)
            poroline=[]

        elif (len(line)==0 or "***" in line):
            break
        else:
            for w in line.split():
              poroline.append(float(w))
    poro.append(poroline)
    poro = np.array(poro)

    return poro

def read_density_water_file(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()

    # === Extract X Coordinates ===
    x_start = lines.index(next(l for l in lines if "R-Direction Node Coordinate" in l)) + 3
    x_coords = []
    while lines[x_start].strip():
        parts = lines[x_start].split()
        x_coords.extend([float(p) for p in parts])
        x_start += 1
    x_start = lines.index(next(l for l in lines if "R-Direction Node Coordinate" in l)) + 8
    while lines[x_start].strip():
        parts = lines[x_start].split()
        x_coords.extend([float(p) for p in parts])
        x_start += 1

    # === Extract Y Coordinates ===
    try :
      y_start = lines.index(next(l for l in lines if "Y-Direction Node Coordinate" in l)) + 2
      y_coords = []
      while lines[y_start].strip():
        parts = lines[y_start].split()
        y_coords.extend([float(p) for p in parts])
        y_start += 1
    except:
        y_coords = [1]

    # === Extract Z Coordinates ===
    z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 3
    z_coords = []
    while lines[z_start].strip():
        parts = lines[z_start].split()
        z_coords.extend([float(p) for p in parts])
        z_start += 1
    z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 8
    while lines[z_start].strip():
        parts = lines[z_start].split()
        z_coords.extend([float(p) for p in parts])
        z_start += 1
    # z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 13
    # while lines[z_start].strip():
    #     parts = lines[z_start].split()
    #     z_coords.extend([float(p) for p in parts])
    #     z_start += 1
    # z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 18
    # while lines[z_start].strip():
    #     parts = lines[z_start].split()
    #     z_coords.extend([float(p) for p in parts])
    #     z_start += 1

    # === Extract Temperature Field ===
    temp_start =find_last_occurrence(filepath, "Water Density")+3
    temperature = [] 
    for line in lines[temp_start:]:
        parts = line.strip().split()
        if len(parts) < 2:
            break
        row = [float(p) for p in parts[1:]]  # skip Z-index!exitr
        temperature.append(row)
    for i in range(len(temperature)):
        line = lines[temp_start+4+len(z_coords)+i]
        parts = line.strip().split()
        if len(parts) < 2:
            break
        temperature[i]=np.array(temperature[i]+[float(p) for p in parts[1:]])  # skip Z-index))
        # temperature.append(row)

    # Convert everything to numpy arrays
    x = np.array(x_coords)
    y = np.array(y_coords)
    z = np.array(z_coords)
    T = np.array(temperature)  # shape should be (len(z), len(x))

    return x, y, z, T 

def read_density_steam_file(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()

    # === Extract X Coordinates ===
    x_start = lines.index(next(l for l in lines if "R-Direction Node Coordinate" in l)) + 3
    x_coords = []
    while lines[x_start].strip():
        parts = lines[x_start].split()
        x_coords.extend([float(p) for p in parts])
        x_start += 1
    x_start = lines.index(next(l for l in lines if "R-Direction Node Coordinate" in l)) + 8
    while lines[x_start].strip():
        parts = lines[x_start].split()
        x_coords.extend([float(p) for p in parts])
        x_start += 1

    # === Extract Y Coordinates ===
    try:
      y_start = lines.index(next(l for l in lines if "Y-Direction Node Coordinate" in l)) + 2
      y_coords = []
      while lines[y_start].strip():
        parts = lines[y_start].split()
        y_coords.extend([float(p) for p in parts])
        y_start += 1
    except:
        y_coords = [1]

    # === Extract Z Coordinates ===
    z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 3
    z_coords = []
    while lines[z_start].strip():
        parts = lines[z_start].split()
        z_coords.extend([float(p) for p in parts])
        z_start += 1
    z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 8
    while lines[z_start].strip():
        parts = lines[z_start].split()
        z_coords.extend([float(p) for p in parts])
        z_start += 1
    # z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 13
    # while lines[z_start].strip():
    #     parts = lines[z_start].split()
    #     z_coords.extend([float(p) for p in parts])
    #     z_start += 1
    # z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 18
    # while lines[z_start].strip():
    #     parts = lines[z_start].split()
    #     z_coords.extend([float(p) for p in parts])
    #     z_start += 1

    # === Extract Temperature Field ===
    temp_start =find_last_occurrence(filepath, "Steam Density")+3
    temperature = [] 
    for line in lines[temp_start:]:
        parts = line.strip().split()
        if len(parts) < 2:
            break
        row = [float(p) for p in parts[1:]]  # skip Z-index!exitr
        temperature.append(row)
    for i in range(len(temperature)):
        line = lines[temp_start+4+len(z_coords)+i]
        parts = line.strip().split()
        if len(parts) < 2:
            break
        temperature[i]=np.array(temperature[i]+[float(p) for p in parts[1:]])  # skip Z-index))
        # temperature.append(row)

    # Convert everything to numpy arrays
    x = np.array(x_coords)
    y = np.array(y_coords)
    z = np.array(z_coords)
    T = np.array(temperature)  # shape should be (len(z), len(x))

    return x, y, z, T 


def read_pressure_file(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()

    # === Extract X Coordinates ===
    x_start = lines.index(next(l for l in lines if "R-Direction Node Coordinate" in l)) + 3
    x_coords = []
    while lines[x_start].strip():
        parts = lines[x_start].split()
        x_coords.extend([float(p) for p in parts])
        x_start += 1
    x_start = lines.index(next(l for l in lines if "R-Direction Node Coordinate" in l)) + 8
    while lines[x_start].strip():
        parts = lines[x_start].split()
        x_coords.extend([float(p) for p in parts])
        x_start += 1

    # === Extract Y Coordinates ===
    try:
      y_start = lines.index(next(l for l in lines if "Y-Direction Node Coordinate" in l)) + 2
      y_coords = []
      while lines[y_start].strip():
        parts = lines[y_start].split()
        y_coords.extend([float(p) for p in parts])
        y_start += 1
    except:
        y_coords = [1]

    # === Extract Z Coordinates ===
    z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 3
    z_coords = []
    while lines[z_start].strip():
        parts = lines[z_start].split()
        z_coords.extend([float(p) for p in parts])
        z_start += 1
    z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 8
    while lines[z_start].strip():
        parts = lines[z_start].split()
        z_coords.extend([float(p) for p in parts])
        z_start += 1
    # z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 13
    # while lines[z_start].strip():
    #     parts = lines[z_start].split()
    #     z_coords.extend([float(p) for p in parts])
    #     z_start += 1
    # z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 18
    # while lines[z_start].strip():
    #     parts = lines[z_start].split()
    #     z_coords.extend([float(p) for p in parts])
    #     z_start += 1

    # === Extract Temperature Field ===
    temp_start =find_last_occurrence(filepath, "Pressure")+3
    temperature = [] 
    for line in lines[temp_start:]:
        parts = line.strip().split()
        if len(parts) < 2:
            break
        row = [float(p) for p in parts[1:]]  # skip Z-index!exitr
        temperature.append(row)
    for i in range(len(temperature)):
       
        line = lines[temp_start+4+len(z_coords)+i]
        parts = line.strip().split()
        if len(parts) < 2:
            break
        temperature[i]=np.array(temperature[i]+[float(p) for p in parts[1:]])  # skip Z-index))
        # temperature.append(row)

    # Convert everything to numpy arrays
    x = np.array(x_coords)
    y = np.array(y_coords)
    z = np.array(z_coords)
    T = np.array(temperature)  # shape should be (len(z), len(x))

    return x, y, z, T




def read_watersat_file(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()

    # === Extract X Coordinates ===
    x_start = lines.index(next(l for l in lines if "R-Direction Node Coordinate" in l)) + 3
    x_coords = []
    while lines[x_start].strip():
        parts = lines[x_start].split()
        x_coords.extend([float(p) for p in parts])
        x_start += 1
    x_start = lines.index(next(l for l in lines if "R-Direction Node Coordinate" in l)) + 8
    while lines[x_start].strip():
        parts = lines[x_start].split()
        x_coords.extend([float(p) for p in parts])
        x_start += 1

    # === Extract Y Coordinates ===
    try:
      y_start = lines.index(next(l for l in lines if "Y-Direction Node Coordinate" in l)) + 2
      y_coords = []
      while lines[y_start].strip():
        parts = lines[y_start].split()
        y_coords.extend([float(p) for p in parts])
        y_start += 1
    except:
        y_coords = [1]

    # === Extract Z Coordinates ===
    z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 3
    z_coords = []
    while lines[z_start].strip():
        parts = lines[z_start].split()
        z_coords.extend([float(p) for p in parts])
        z_start += 1
    z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 8
    while lines[z_start].strip():
        parts = lines[z_start].split()
        z_coords.extend([float(p) for p in parts])
        z_start += 1
    # z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 13
    # while lines[z_start].strip():
    #     parts = lines[z_start].split()
    #     z_coords.extend([float(p) for p in parts])
    #     z_start += 1
    # z_start = lines.index(next(l for l in lines if "Z-Direction Node Coordinate" in l)) + 18
    # while lines[z_start].strip():
    #     parts = lines[z_start].split()
    #     z_coords.extend([float(p) for p in parts])
    #     z_start += 1

    # === Extract Saturation Field ===
    temp_start =find_last_occurrence(filepath, "Water Saturation")+3
    saturation = [] 
    for line in lines[temp_start:]:
        parts = line.strip().split()
        if len(parts) < 2:
            break
        row = [float(p) for p in parts[1:]]  # skip Z-index!exitr
        saturation.append(row)
    for i in range(len(saturation)):
        line = lines[temp_start+4+len(z_coords)+i]
        parts = line.strip().split()
        if len(parts) < 2:
            break
        saturation[i]=np.array(saturation[i]+[float(p) for p in parts[1:]])  # skip Z-index))
        # temperature.append(row)
    
    # === Extract mass fraction Field ===
    temp_start =find_last_occurrence(filepath, "Mass Fraction")+3
    massfra = [] 
    for line in lines[temp_start:]:
        parts = line.strip().split()
        if len(parts) < 2:
            break
        row = [float(p) for p in parts[1:]]  # skip Z-index!exitr
        massfra.append(row)
    for i in range(len(massfra)):
        line = lines[temp_start+4+len(z_coords)+i]
        parts = line.strip().split()
        if len(parts) < 2:
            break
        massfra[i]=np.array(massfra[i]+[float(p) for p in parts[1:]])  # skip Z-index))
        # temperature.append(row)

    # Convert everything to numpy arrays
    massfra = np.array(massfra)  # shape should be (len(z), len(x))

    # Convert everything to numpy arrays
    x = np.array(x_coords)
    y = np.array(y_coords)
    z = np.array(z_coords)
    satu = np.array(saturation)  # shape should be (len(z), len(x))

    return x, y, z, saturation, massfra




def boom_eval(probname):
    x, y, z, temperature = read_temperature_file("Out_temperature."+probname+'.out')
    x, y, z, pressure = read_pressure_file("Out_pressure."+probname+'.out')
    pressure *=1e-4 #from dyne/cm2 to kPa
    # print("qui A abbiamo x e z:", len(x), len(z), temperature.shape)
    x, y, z, waterrho = read_density_water_file("Out_density."+probname+'.out')
    x, y, z, watervolsat, watermasfra = read_watersat_file("Out_saturation."+probname+'.out')
    watervolsat=np.array(watervolsat)
    Sl = watervolsat
    Sg = 1-watervolsat
    
    # print("qui B abbiamo x e z:", len(x), len(z), waterrho.shape)
    x, y, z, steamrho = read_density_steam_file("Out_density."+probname+'.out')
    poro = read_porosity("Out_porosity."+probname+'.out')
    # print("qui C abbiamo x e z:", len(x), len(z))
    rockrho = read_rock_density(probname+'.in')
    endtime = read_sim_end(probname)



    temperature = temperature + 273.15
    waterrho *= 1000 
    steamrho *= 1000

    S3 = 2180   #kPa
    Pthresh = 2*S3*(1-poro)/(3*poro*np.sqrt(poro**(-1/3)-1))
    
    overcome = pressure/Pthresh
    overcome_count = np.sum(overcome[1:,:]>1)/len(np.ravel(overcome[1:,:]))



    #boom_total = eos.boom_rev(steamrho, Sg, waterrho, Sl, temperature)
    print("Calculating monophase explosivity...")
    boom_total_1p = eos.boom_irr_monophase(steamrho, waterrho, Sg, Sl, temperature)
    print("Calculating biphase explosivity")
    boom_total_2p = eos.boom_irr_biphase  (steamrho, waterrho, Sg, Sl, temperature, pressure)
    print("Calculating final explosivity...")
    #total weight is the average between monophase and biphase, 
    #but if biphase  is nan then we set to 0
    boom_2pweight = 1-(np.isnan(boom_total_2p))  # wheight is 1 if not nan
    boom_1pweight = 1-(np.isnan(boom_total_1p))  #
    #do average
    boom_total = (boom_2pweight*np.nan_to_num(boom_total_2p) +  boom_1pweight*np.nan_to_num(boom_total_1p))/(boom_2pweight+boom_1pweight)
    integrated_explo = integrate_explo(x, z, boom_total)
    for i in range(len(boom_total[0,:])):
        boom_total[0,i] = 0
        
    # rock_potential_energy = integrate_rock_weight(x,z,rockrho)
    # print("qui D abbiamo x e z:", len(x), len(z))
    # print("Explosive energy:", integrated_explo)
    # print("Rock potential energy:", rock_potential_energy)
    #print("integrated explo",integrated_explo)
    #print("rock_potential_energy", rock_potential_energy)

    return overcome_count, endtime, integrated_explo


def integrate_rock_weight(x,z,rockdensi):
    return (x[-1]-x[0])*(z[-1])*rockdensi*9.81*1000
    
def integrate_explo(x, z, boom_total):
    nx = len(x)
    nz = len(z)
    boom_integrated = 0
    #print(boom_total.shape, nx, nz)
    #print(x)
    #print(z)
    for i in range(nx-1):
        for j in range(nz-1):
            boom_integrated += (x[i+1]-[i])*(z[j+1]-z[j])*boom_total[j,i]
    return boom_integrated


def read_rock_density(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()
    iline = 0
    for line in lines:
        if ("rxden" in line):
            iline+=2
            break
        else:
            iline+=1
    den_string = lines[iline].split()[1]
    if ("," in den_string):
        den_string = den_string[:-1]
    rockdensity = float(den_string)
    return rockdensity 

def read_sim_end(probname):
    filepath = "Calc_log."+probname+'.out'
    with open(filepath, 'r') as fin:
        lines = fin.readlines()
    for line in lines[::-1]:
        if ("Last time value calculated") in line:
            wrds = line.split()
            endtime = float(wrds[5])
            break
        elif ("Simulation time" ) in line:
            wrds = line.split()
            endtime = float(wrds[-2])
            break
    return endtime

# import os
# #def integrate_weight
# os.chdir("tuno0000")
# print(boom_eval("tuno0000"))

    
