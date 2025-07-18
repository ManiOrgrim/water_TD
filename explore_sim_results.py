import numpy as np

mago = np.load("tqat.npy")


import matplotlib.pyplot as plt

fag, ax = plt.subplots(dpi=300) 
ax.scatter(mago[:,1], mago[:,-1], s =4)
ax.set_title("Flux vs explosive energy")
ax.set_xlabel(r"Flux [J/m$^2$ d]")
ax.set_ylabel(r"Explosive energy [J]")

fbg, bx = plt.subplots(dpi=300) 
bx.scatter(mago[:,1], mago[:,-2], s =4)
bx.set_title("Flux vs endtime")
bx.set_xlabel(r"Flux [J/m$^2$ d]")
bx.set_ylabel(r"Endtime [years]")

fcg, cx = plt.subplots(dpi=300) 
cx.scatter(mago[:,1], mago[:,-3], s =4)
cx.set_title("Flux vs MOP")
cx.set_xlabel(r"Flux [J/m$^2$ d]")
cx.set_ylabel(r"MOP [-]")
cx.set_yscale("log")


fdg, dx = plt.subplots(dpi=300) 
dx.scatter(mago[:,2], mago[:,-1], s =4)
dx.set_title(r"T$_b$ vs explosive energy")
dx.set_xlabel(r"T$_b$ [K]")
dx.set_ylabel(r"Explosive energy [J]")

feg, ex = plt.subplots(dpi=300) 
ex.scatter(mago[:,2], mago[:,-2], s =4)
ex.set_title(r"T$_b$vs endtime")
ex.set_xlabel(r"T$_b$ [K]")
ex.set_ylabel(r"Endtime [years]")

ffg, fx = plt.subplots(dpi=300) 
fx.scatter(mago[:,2], mago[:,-3], s =4)
fx.set_title(r"T$_b$ vs MOC")
fx.set_xlabel(r"T$_b$ [K]")
fx.set_ylabel(r"MOP [-]")
fx.set_yscale("log")


fdg, dx = plt.subplots(dpi=300) 
dx.scatter(mago[:,3], mago[:,-1], s =4)
dx.set_title(r"High permeability vs explosive energy")
dx.set_xlabel(r"High permeability [m$^2$]")
dx.set_ylabel(r"Explosive energy [J]")

feg, ex = plt.subplots(dpi=300) 
ex.scatter(mago[:,3], mago[:,-2], s =4)
ex.set_title(r"High permeability vs endtime")
ex.set_xlabel(r"High permeability [m$^2$]")
ex.set_ylabel(r"Endtime [years]")

fgg, gx = plt.subplots(dpi=300) 
gx.scatter(mago[:,3], mago[:,-3], s =4)
gx.set_title(r"High permeability vs MOP")
gx.set_xlabel(r"High permeability [m$^2$]")
gx.set_ylabel(r"MOP [-]")
gx.set_yscale("log")


fhg, hx = plt.subplots(dpi=300) 
hx.scatter(mago[:,4], mago[:,-1], s =4)
hx.set_title(r"High porosity vs explosive energy")
hx.set_xlabel(r"High porosity [-]")
hx.set_ylabel(r"Explosive energy [J]")

fig, ix = plt.subplots(dpi=300) 
ix.scatter(mago[:,4], mago[:,-2], s =4)
ix.set_title(r"High porosity vs endtime")
ix.set_xlabel(r"High porosity [-]")
ix.set_ylabel(r"Endtime [years]")

fjg, jx = plt.subplots(dpi=300) 
jx.scatter(mago[:,4], mago[:,-3], s =4)
jx.set_title(r"High porosity vs MOP")
jx.set_xlabel(r"High porosity [-]")
jx.set_ylabel(r"MOP [-]")
jx.set_yscale("log")



fkg, kx = plt.subplots(dpi=300) 
kx.scatter(mago[:,5], mago[:,-1], s =4)
kx.set_title(r"Low permeability vs explosive energy")
kx.set_xlabel(r"Low  permeability [m$^2$]")
kx.set_ylabel(r"Explosive energy [J]")

flg, lx = plt.subplots(dpi=300) 
lx.scatter(mago[:,5], mago[:,-2], s =4)
lx.set_title(r"Low  permeability vs endtime")
lx.set_xlabel(r"Low  permeability [m$^2$]")
lx.set_ylabel(r"Endtime [years]")

fmg, mx = plt.subplots(dpi=300) 
mx.scatter(mago[:,5], mago[:,-3], s =4)
mx.set_title(r"Low  permeability vs MOP")
mx.set_xlabel(r"Low  permeability [m$^2$]")
mx.set_ylabel(r"MOP [-]")
mx.set_yscale("log")


fng, nx = plt.subplots(dpi=300) 
nx.scatter(mago[:,6], mago[:,-1], s =4)
nx.set_title(r"Low  porosity vs explosive energy")
nx.set_xlabel(r"Low  porosity [-]")
nx.set_ylabel(r"Explosive energy [J]")

fog, ox = plt.subplots(dpi=300) 
ox.scatter(mago[:,6], mago[:,-2], s =4)
ox.set_title(r"Low  porosity vs endtime")
ox.set_xlabel(r"Low  porosity [-]")
ox.set_ylabel(r"Endtime [years]")

fpg, px = plt.subplots(dpi=300) 
px.scatter(mago[:,6], mago[:,-3], s =4)
px.set_title(r"Low  porosity vs MOP")
px.set_xlabel(r"Low  porosity [-]")
px.set_ylabel(r"MOP [-]")
px.set_yscale("log")



fqg, qx = plt.subplots(dpi=300) 
qx.scatter(mago[:,7], mago[:,-1], s =4)
qx.set_title(r"Dz vs explosive energy")
qx.set_xlabel(r"Dz[m]")
qx.set_ylabel(r"Explosive energy [J]")

frg, rx = plt.subplots(dpi=300) 
rx.scatter(mago[:,7], mago[:,-2], s =4)
rx.set_title(r"Dz vs endtime")
rx.set_xlabel(r"Dz [m]")
rx.set_ylabel(r"Endtime [years]")

fsg, sx = plt.subplots(dpi=300) 
sx.scatter(mago[:,7], mago[:,-3], s =4)
sx.set_title(r"Dz vs MOP")
sx.set_xlabel(r"Dz [m]")
sx.set_ylabel(r"MOP [-]")
sx.set_yscale("log")




ftg, tx = plt.subplots(dpi=300) 
tx.scatter(mago[:,-1], mago[:,-2], s =4)
tx.set_title(r"Explosive energy vs endtime")
tx.set_xlabel(r"Expolosive energy [J]")
tx.set_ylabel(r"Endtime [years]")

fug, ux = plt.subplots(dpi=300) 
ux.scatter(mago[:,-1], mago[:,-3], s =4)
ux.set_title(r"Explosive energy vs MOP")
ux.set_xlabel(r"Explosive energy [J]")
ux.set_ylabel(r"MOP [-]")
ux.set_yscale("log")
ux.set_xscale("log")


fvg, vx = plt.subplots(dpi=300) 
vx.scatter(mago[:,-2], mago[:,-3], s =4)
vx.set_title(r"Endtime vs MOP")
vx.set_xlabel(r"Endtime [years]")
vx.set_ylabel(r"MOP [-]")
vx.set_yscale("log")

