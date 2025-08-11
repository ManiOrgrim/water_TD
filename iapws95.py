# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 07:54:20 2025

@author: manim
"""

import numpy as np

Tc = 647.096
rho_c = 322.0          
R = 0.46151805       

Ni0 =np.array( [
    -8.32044648201,
     6.68321052680,
     3.00632000000,
     0.01243600000,
     0.97315000000,
     1.27950000000,
     0.96956000000,
     0.24873000000
])

Gammai0 = np.array([
    0.0,
    0.0,
    0.0,
    1.28728967,
    3.53734222,
    7.74073708,
    9.24437796,
    27.50751050
])

Ni = np.array([
    0.12533547935523e-1,
    0.78957634722828e1,
   -0.87803203303561e1,
    0.31802509345418,
   -0.26145533859358,   #5
   -0.78199751687981e-2,
    0.88089493102134e-2,
   -0.66856572307965,
    0.20433810950965,
   -0.66212605039687e-4, #10
   -0.19232721156002, 
   -0.25709043003438,
    0.16074868486251,
   -0.40092828925807e-1,
    0.39343422603254e-6, #15
   -0.75941377088144e-5,
    0.56250979351888e-3,
   -0.15608652257135e-4,
    0.11537996422951e-8,
    0.36582165144204e-6, #20
   -0.13251180074668e-11, 
   -0.62639586912454e-9,
   -0.10793600908932,
    0.17611491008752e-1,
    0.22132295167546,#25
   -0.40247669763528,
    0.58083399985759,
    0.49969146990806e-2,
   -0.31358700712549e-1,
   -0.74315929710341, #30
    0.47807329915480,
    0.20527940895948e-1,
   -0.13636435110343,
    0.14180634400617e-1,
    0.83326504880713e-2,  #35
   -0.29052336009585e-1,
    0.38615085574206e-1,
   -0.20393486513704e-1,
   -0.16554050063734e-2,
    0.19955571979541e-2, #40
    0.15870308324157e-3,
   -0.16388568342530e-4,
    0.43613615723811e-1,
    0.34994005463765e-1,
   -0.76788197844621e-1,  #45
    0.22446277332006e-1,
   -0.62689710414685e-4,
   -0.55711118565645e-9,
   -0.19905718354408,
    0.31777497330738, #50
   -0.11841182425981,
   -0.31306260323435e2,
    0.31546140237781e2,
   -0.25213154341695e4,
   -0.14874640856724,  #55
    0.31806110878444
])

Ci = np.array([   # GIUSTO
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  # 1–7
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  # 8-14 
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  # 15-21
    1.0,  # 22
    2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,  # 23-29
    2.0, 2.0, 2.0,  #30-32
    2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,  # 32–42
    3.0, 3.0, 3.0, 3.0,  # 43–46
    4.0,                # 47
    6.0, 6.0, 6.0, 6.0,  # 48–51
    0.0, 0.0, 0.0, 0.0, 0.0  # 52–56 (restano a 0.0)
])

Di = np.array([   #GIUSTO
    1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 4.0, 1.0, 1.0, 1.0,
    2.0, 2.0, 3.0, 4.0, 4.0, 5.0, 7.0, 9.0, 10.0, 11.0,
    13.0, 15.0, 1.0, 2.0, 2.0, 2.0, 3.0, 4.0, 4.0, 4.0,
    5.0, 6.0, 6.0, 7.0, 9.0, 9.0, 9.0, 9.0, 9.0, 10.0,
    10.0, 12.0, 3.0, 4.0, 4.0, 5.0, 14.0, 3.0, 6.0, 6.0,
    6.0, 3.0, 3.0, 3.0, 0.0, 0.0  # gli ultimi 2 (55, 56) restano a 0.0
])

Ti = np.array([ #GIUSTO
   -0.5, 0.875, 1.0, 0.5, 0.75, 0.375, 1.0,
    4.0, 6.0, 12.0, 1.0, 5.0, 4.0, 2.0, 13.0,
    9.0, 3.0, 4.0, 11.0, 4.0, 13.0, 1.0, 7.0,
    1.0, 9.0, 10.0, 10.0, 3.0, 7.0, 10.0, 10.0,
    6.0, 10.0, 10.0, 1.0, 2.0, 3.0, 4.0, 8.0,
    6.0, 9.0, 8.0, 16.0, 22.0, 23.0, 23.0, 10.0,
    50.0, 44.0, 46.0, 50.0, 0.0, 1.0, 4.0, 0.0
])

alphai = [0.0] * 56
alphai[51] = 20.0
alphai[52] = 20.0
alphai[53] = 20.0
alphai=np.array(alphai)

betai = [0.0]*56
betai = np.array(betai)
betai[51] = 150
betai[52] = 150
betai[53] = 250

gammai = [0.0]*56
gammai[51] = 1.21
gammai[52] = 1.21
gammai[53] = 1.25
gammai = np.array(gammai)

Epsi = [0.0]*56
Epsi[51] = 1
Epsi[52] = 1
Epsi[53] = 1
Epsi= np.array(Epsi)

def get_dPhi0_tau(delta, tau):
    dphi0_tau = Ni0[1]+Ni0[2]/tau
    for i in range(3,7):
        dphi0_tau = dphi0_tau+Ni0[i]*Gammai0[i]*((1.0-np.exp(-Gammai0[i]*tau))**(-1.0)-1.0)
    return  dphi0_tau 

def get_dPhiR_tau(delta, tau):
    # Assumiamo che queste funzioni restituiscano array numpy di lunghezza 56
    dPhiR_tau = 0.0

    # Primo ciclo: i = 0 a 6 (Fortran 1–7 → Python 0–6)
    for i in range(0, 7):
        dPhiR_tau += Ni[i] * Ti[i] * (delta ** Di[i]) * tau ** (Ti[i] - 1.0)

    # Secondo ciclo: i = 7 a 50 (Fortran 8–51 → Python 7–50)
    for i in range(7, 51):
        dPhiR_tau += Ni[i] * Ti[i] * (delta ** Di[i]) * tau ** (Ti[i] - 1.0) * np.exp(-delta ** Ci[i])

    # Terzo ciclo: i = 51 a 53 (Fortran 52–54 → Python 51–53)
    for i in range(51, 54):
        exp_term = np.exp(-alphai[i] * (delta - Epsi[i])**2 - betai[i] * (tau - gammai[i])**2)
        part = Ti[i] / tau - 2.0 * betai[i] * (tau - gammai[i])
        dPhiR_tau += Ni[i] * (delta ** Di[i]) * (tau ** Ti[i]) * exp_term * part

    # Primo gruppo di costanti per il termine 55 (indice Python 54)
    Aii = 0.32
    betaii = 0.3
    Bii = 0.2
    apiccoloii = 3.5
    Cii = 28.0
    Dii = 700.0
    bpiccoloii = 0.85

    theta = (1.0 - tau) + Aii * ((delta - 1.0) ** 2) ** (1.0 / (2.0 * betaii))
    Delta_cap = theta**2 + Bii * ((delta - 1.0) ** 2) ** apiccoloii
    argexp = -Cii * (delta - 1.0)**2 - Dii * (tau - 1.0)**2
    Psi = np.exp(argexp)

    dPsi_tau = -2.0 * Dii * (tau - 1.0) * Psi
    dDeltacapbi_tau = -2.0 * theta * bpiccoloii * (Delta_cap ** (bpiccoloii - 1.0))
    dPhiR_tau += Ni[54] * delta * (dDeltacapbi_tau * Psi + (Delta_cap ** bpiccoloii) * dPsi_tau)

    # Secondo gruppo di costanti per il termine 56 (indice Python 55)
    Aii = 0.32
    betaii = 0.3
    Bii = 0.2
    apiccoloii = 3.5
    Cii = 32.0
    Dii = 800.0
    bpiccoloii = 0.95

    theta = (1.0 - tau) + Aii * ((delta - 1.0) ** 2) ** (1.0 / (2.0 * betaii))
    Delta_cap = theta**2 + Bii * ((delta - 1.0) ** 2) ** apiccoloii
    argexp = -Cii * (delta - 1.0)**2 - Dii * (tau - 1.0)**2
    Psi = np.exp(argexp)

    dPsi_tau = -2.0 * Dii * (tau - 1.0) * Psi
    dDeltacapbi_tau = -2.0 * theta * bpiccoloii * (Delta_cap ** (bpiccoloii - 1.0))
    dPhiR_tau += Ni[55] * delta * (dDeltacapbi_tau * Psi + (Delta_cap ** bpiccoloii) * dPsi_tau)

    return dPhiR_tau

def get_dPhiR_delta(delta, tau):
    dPhir_delta = 0.0

    # Prima parte: i=1..7 (Fortran 1-based)
    for i in range(7):  # python 0-based indices 0..6
        dPhir_delta += Ni[i] * Di[i] * delta**(Di[i] - 1.0) * tau**Ti[i]
 #       print("L1", dPhir_delta)
#    print("L1", dPhir_delta)

    # Parte i=8..51 (python 7..50)
    for i in range(7, 51):
        argexp = delta**Ci[i]

        dPhir_delta += Ni[i] * (delta**(Di[i] - 1.0)) * (tau**Ti[i]) * \
                       (Di[i] - Ci[i] * delta**Ci[i]) * np.exp(-argexp)

    # Parte i=52..54 (python 51..53)
    for i in range(51, 54):
        exp_arg = -alphai[i] * (delta - Epsi[i])**2 - betai[i] * (tau - gammai[i])**2
        dPhir_delta += Ni[i] * (delta**Di[i]) * (tau**Ti[i]) * np.exp(exp_arg) * \
                       (Di[i] / delta - 2.0 * alphai[i] * (delta - Epsi[i]))

    # Costanti
    Aii = 0.32
    betaii = 0.3
    Bii = 0.2
    apiccoloii = 3.5
    Cii = 28.0
    Dii = 700.0
    bpiccoloii = 0.85

    theta = (1.0 - tau) + Aii * ((delta - 1.0)**2)**(1.0 / (2.0 * betaii))
    Delta_cap = theta**2 + Bii * ((delta - 1.0)**2)**apiccoloii
    argexp = -Cii * (delta - 1.0)**2 - Dii * (tau - 1.0)**2

    Psi = np.exp(argexp)

    dDeltacap_delta = (delta - 1.0) * (Aii * theta * (2.0 / betaii) * \
                        ((delta - 1.0)**2)**(1.0 / (2.0 * betaii) - 1.0) + \
                        2.0 * Bii * apiccoloii * ((delta - 1.0)**2)**(apiccoloii - 1.0))

    dPsi_delta = -2.0 * Cii * (delta - 1.0) * Psi
    dDeltacapbi_delta = bpiccoloii * Delta_cap**(bpiccoloii - 1.0) * dDeltacap_delta

    dPhir_delta += Ni[54] * ((Psi + delta * dPsi_delta) * Delta_cap**bpiccoloii + \
                             dDeltacapbi_delta * delta * Psi)

    # Seconda parte con costanti diverse
    Cii = 32.0
    Dii = 800.0
    bpiccoloii = 0.95

    theta = (1.0 - tau) + Aii * ((delta - 1.0)**2)**(1.0 / (2.0 * betaii))
    Delta_cap = theta**2 + Bii * ((delta - 1.0)**2)**apiccoloii
    argexp = -Cii * (delta - 1.0)**2 - Dii * (tau - 1.0)**2

    Psi = np.exp(argexp)

    dDeltacap_delta = (delta - 1.0) * (Aii * theta * (2.0 / betaii) * \
                        ((delta - 1.0)**2)**(1.0 / (2.0 * betaii) - 1.0) + \
                        2.0 * Bii * apiccoloii * ((delta - 1.0)**2)**(apiccoloii - 1.0))

    dPsi_delta = -2.0 * Cii * (delta - 1.0) * Psi
    dDeltacapbi_delta = bpiccoloii * Delta_cap**(bpiccoloii - 1.0) * dDeltacap_delta

    dPhir_delta += Ni[55] * ((Psi + delta * dPsi_delta) * Delta_cap**bpiccoloii + \
                             dDeltacapbi_delta * delta * Psi)

    return dPhir_delta

def get_phi_ideal(delta, tau):
    # Calcola phi0 (funzione ideale)
    phi0 = np.log(delta) + Ni0[0] + Ni0[1]*tau + Ni0[2]*np.log(tau)

    for i in range(3, 8):  # in Python gli indici partono da 0
        phi0 += Ni0[i] * np.log(1.0 - np.exp(-Gammai0[i] * tau))
    
    return phi0


def get_phi_residual(delta, tau):
    # Inizializzo gli array
    
    phir = 0.0

    # Primo ciclo 1-7 (Python 0-6)
    for i in range(7):
        phir += Ni[i] * (delta ** Di[i]) * (tau ** Ti[i])

    # Secondo ciclo 8-51 (Python 7-50)
    for i in range(7, 51):
        exparg = - (delta ** Ci[i])
        phir += Ni[i] * (delta ** Di[i]) * (tau ** Ti[i]) * np.exp(exparg)

    # Terzo ciclo 52-54 (Python 51-53)
    for i in range(51, 54):
        phir += Ni[i] * (delta ** Di[i]) * (tau ** Ti[i]) * np.exp(
            -alphai[i] * (delta - Epsi[i])**2 - betai[i] * (tau - gammai[i])**2
        )

    # Costanti
    Aii = 0.32
    betaii = 0.3
    Bii = 0.2
    apiccoloii = 3.5
    Cii = 28.0
    Dii = 700.0
    bpiccoloii = 0.85

    theta = (1.0 - tau) + Aii * ((delta - 1.0)**2) ** (1.0 / (2.0 * betaii))
    Delta_cap = theta**2 + Bii * ((delta - 1.0)**2) ** apiccoloii
    exparg = -Cii * (delta - 1.0)**2 - Dii * (tau - 1.0)**2
    Psi = np.exp(exparg)

    phir += Ni[54] * (Delta_cap ** bpiccoloii) * delta * Psi
    
    # Secondo gruppo di costanti
    Cii = 32.0
    Dii = 800.0
    bpiccoloii = 0.95

    theta = (1.0 - tau) + Aii * ((delta - 1.0)**2) ** (1.0 / (2.0 * betaii))
    Delta_cap = theta**2 + Bii * ((delta - 1.0)**2) ** apiccoloii
    exparg = -Cii * (delta - 1.0)**2 - Dii * (tau - 1.0)**2
    Psi = np.exp(exparg)

    phir += Ni[55] * (Delta_cap ** bpiccoloii) * delta * Psi

    return phir


def EOS_intenergy(rho, T):
    delta = rho/rho_c
    tau   = Tc/T
    dphi0_tau = get_dPhi0_tau(delta, tau)
    dphiR_tau = get_dPhiR_tau(delta,tau)
    intenergy = R*T*tau*(dphi0_tau+dphiR_tau)
    return intenergy

def EOS_pressure(rho, T):
    delta = rho/rho_c
    tau = Tc/T
    dPhiR_delta = get_dPhiR_delta(delta, tau)
    P = rho*R*T*(1.0+delta*dPhiR_delta)
    return P

def EOS_entropy(rho,T):
    delta = rho/rho_c
    tau = Tc/T
    phi0 = get_phi_ideal(delta, tau)
    phiR = get_phi_residual(delta, tau)
    dPhi0_tau = get_dPhi0_tau(delta, tau)
    dPhiR_tau = get_dPhiR_tau(delta, tau)
    S = R*(tau*(dPhi0_tau+dPhiR_tau)-phi0-phiR)
    return S


def EOS_gibbs(rho, T):
    delta = rho / rho_c
    tau = Tc / T
    phi0 = get_phi_ideal(delta, tau)
    phiR = get_phi_residual(delta, tau)
    dPhiR_delta = get_dPhiR_delta(delta, tau)

    Gibbs = R * T * (1.0 + phi0 + phiR + delta * dPhiR_delta)
    return Gibbs


def find_equilibrium_broyden(T):
    """
    Trova l'equilibrio di fase tramite il metodo di Broyden.
    
    Parametri:
    T           - temperatura (float)
    rho_g_init  - densità gas iniziale (float, input/output)
    rho_l_init  - densità liquido iniziale (float, input/output)
    Tc          - temperatura critica (float)
    
    Restituisce:
    rho_g, rho_l, P
    """
    # Inizializza xprev in base a T
    #print(T)
    if T > Tc:
        raise ValueError("T maggiore della temperatura critica: caso non gestito.")
    elif T > 645.0:
        xprev = np.array([240.0, 400.0])
    elif T > 640.0:
        xprev = np.array([100.0, 600.0])
    elif T > 600.0:
        xprev = np.array([50.0, 700.0])
    elif T > 455.0:
        xprev = np.array([10.0, 900.0])
    elif T > 363.0:
        xprev = np.array([1.0, 1000.0])
    else:
        xprev = np.array([0.01, 1000.0])

    epsil = 1e-6


    # Calcola Jacobiano 2x2
    Jaco0 = np.zeros((2, 2))
    
    Jaco0[0, 0] = (EOS_pressure(xprev[0]+epsil, T) -EOS_pressure(xprev[0]-epsil,T))/(2.0*epsil)
    Jaco0[0,1] = -(EOS_pressure(xprev[1]+epsil, T) -EOS_pressure(xprev [1]-epsil, T))/(2.0*epsil)
    Jaco0[1,0] = (EOS_gibbs(xprev[0]+epsil, T) -EOS_gibbs(xprev[0]-epsil, T))/(2.0*epsil)
    Jaco0[1,1] = -(EOS_gibbs(xprev[1]+epsil, T) -EOS_gibbs(xprev[1]-epsil, T))/(2.0*epsil)


    det0 = Jaco0[0, 0]*Jaco0[1, 1] - Jaco0[1, 0]*Jaco0[0, 1]
    # print(det0)
    if abs(det0) < 1e-12:
        raise ValueError("Jacobian matrix is singular or near-singular")
        
    Ainv = np.zeros((2, 2))
    Fprev = np.zeros(2)
    Fnow = np.zeros(2)


    Ainv[0,0] = Jaco0[1,1]/det0
    Ainv[1,0]  = -Jaco0[1,0]/det0
    Ainv[0,1]  = -Jaco0[0,1]/det0
    Ainv[1,1]  = Jaco0[0,0]/det0
    Fprev[0]   = EOS_pressure(xprev[0],T)-EOS_pressure(xprev[1],T)
    Fprev[1]   = EOS_gibbs(xprev[0],T)-EOS_gibbs(xprev[1],T)

    xnow = xprev - np.matmul(Ainv, Fprev)

    diff = 1e10

    while diff > 1e-3:
        Fnow[0]   = EOS_pressure(xnow[0],T)-EOS_pressure(xnow[1],T)
        Fnow[1]   = EOS_gibbs(xnow[0],T)-EOS_gibbs(xnow[1],T)
        # print("All'inizio abbiamo abbiamo Fnow, Fprev")
        # print(Fnow, Fprev)
        Y = Fnow - Fprev
        S = xnow - xprev
        diff = np.linalg.norm(S)
        # print("In broyden abbiamo S, AInv, Y")
        # print (S, Ainv, Y)

        sotto = np.dot(np.matmul(S,Ainv),  Y)
        # print("In broyden abbiamo sotto")
        # print(sotto)
        # if abs(sotto) < 1e-12:
        #     raise ValueError("Denominator too small in Broyden update")

        # Aggiornamento matrice inversa (Broyden)
        term = np.dot(S-np.matmul(Ainv, Y), S)/sotto
        # print("term:", term) 
        # print(Ainv.shape, term.shape)
        Ainv = Ainv +Ainv*term

        xprev = xnow
        Fprev[0] = Fnow[0]
        Fprev[1] = Fnow[1]
        xnow = xprev -np.matmul(Ainv, Fprev)
        # print("In broyden abbiamo xnow, xprev")
        # print(xnow, xprev)
        # print("In broyden abbiamo Fnow, Fprev")
        # print(Fnow, Fprev)
        #controlla che questo ciclo cicli bene

    rho_g, rho_l = xnow[0], xnow[1]
    P = EOS_pressure(rho_g, T)

    return rho_g, rho_l, P



def find_equi_TfromP(P):
    T = 280
    Ptest = P+100
    while(abs(P-Ptest)>8e-2):
        rhog, rhol, Ptest =find_equilibrium_broyden(T)
        # print(rhog, rhog)
        T = T+0.01
        # print("In find equi abbiamo T, rhog, rhol, Ptest")
        # print(T, rhog, rhol, P, Ptest)
    return T, rhog, rhol
    

def boom_rev_orig(rho0_g, Satu0_g, rho0_l, Satu0_l, T0):
    Patm = 101.325
    U0 = rho0_g*Satu0_g*EOS_intenergy(rho0_g, T0) + rho0_l*Satu0_l*EOS_intenergy(rho0_l, T0)
    Entro0 = rho0_g*Satu0_g *np.minimum(EOS_entropy(rho0_g, T0),1000000000) + rho0_l*Satu0_l*np.minimum(EOS_entropy(rho0_l, T0),10000000)
    #Teq, rhog_eq, rhol_eq  = find_equi_TfromP(Patm)
    Teq = 373.12
    rhog_eq = 0.5973
    rhol_eq = 958.378
    Entro_g_atm = rhog_eq*EOS_entropy(rhog_eq, Teq)
    Entro_l_atm = rhol_eq*EOS_entropy(rhol_eq, Teq)
    f = (Entro0-Entro_l_atm)/(Entro_g_atm-Entro_l_atm)
    print("Entro0", Entro0)
    print("Entro_g_atm", Entro_g_atm)
    print("Entro_l_atm", Entro_l_atm)
    print(f)    

    # print("In boom rev abbiamo f, Teq, rhog_eq, rhol_eq")
    # print(f, Teq, rhog_eq, rhol_eq)
    U_boom = (rhog_eq*f*EOS_intenergy(rhog_eq, Teq)+rhol_eq*(1-f)*EOS_intenergy(rhol_eq, Teq)-U0)
    for i in range(f.shape[0]):
        for j in range (f.shape[1]):
            if (f[i,j]< 0):
                if (i==0):
                    U_boom[i,j] = 0
                else:
                    rhol_eq, Tf = isoS_expansion_solver(T0[i,j], rho0_l[i,j], Patm)
                    #print("iih")
                    #print(rhol_eq.shape, EOS_intenergy(rhol_eq, Tf).shape)
                    U_boom[i,j] = rhol_eq*EOS_intenergy(rhol_eq, Tf) -U0[i,j] 
            elif (f[i,j]>1):
                rhog_eq, Tf = isoS_expansion_solver(T0[i,j], rho0_g[i,j], Patm)
                U_boom[i,j] = rhog_eq*EOS_intenergy(rhog_eq, Tf) -U0[i,j]

    #print("U_boom", U_boom)
    return U_boom*1000  #EOS calcuates in kJ


def boom_rev(rho0_g, Satu0_g, rho0_l, Satu0_l, T0):
    Patm = 101.325
    U0 = rho0_g*Satu0_g*EOS_intenergy(rho0_g, T0) + rho0_l*Satu0_l*EOS_intenergy(rho0_l, T0)
    Entro0 = rho0_g*Satu0_g *np.minimum(EOS_entropy(rho0_g, T0),1000000000) + rho0_l*Satu0_l*np.minimum(EOS_entropy(rho0_l, T0),10000000)
    #Teq, rhog_eq, rhol_eq  = find_equi_TfromP(Patm)
    Teq = 373.12
    rhog_eq = 0.5973
    rhol_eq = 958.378
    Entro_g_atm = rhog_eq*EOS_entropy(rhog_eq, Teq)
    Entro_l_atm = rhol_eq*EOS_entropy(rhol_eq, Teq)
    
    M0 = rho0_g*Satu0_g + rho0_l*Satu0_l
    #Entro0 già calcolata sopra
    R0 = M0/Entro0
    Sgf = (-rhol_eq*Entro_l_atm*R0 +rhol_eq)/((rhog_eq*Entro_g_atm-rhol_eq*Entro_l_atm)*R0 +rhol_eq-rhog_eq)
        

    # print("In boom rev abbiamo f, Teq, rhog_eq, rhol_eq")
    # print(f, Teq, rhog_eq, rhol_eq)
    U_boom = (rhog_eq*Sgf*EOS_intenergy(rhog_eq, Teq)+rhol_eq*(1-Sgf)*EOS_intenergy(rhol_eq, Teq)-U0)
    for i in range(Sgf.shape[0]):
        for j in range (Sgf.shape[1]):
            if (Sgf[i,j]< 0):
                if (i==0):
                    U_boom[i,j] = 0
                else:
                    rhol_eq, Tf = isoS_expansion_solver(T0[i,j], rho0_l[i,j], Patm)
                    #print("iih")
                    #print(rhol_eq.shape, EOS_intenergy(rhol_eq, Tf).shape)
                    U_boom[i,j] = rhol_eq*EOS_intenergy(rhol_eq, Tf) -U0[i,j] 
            elif (Sgf[i,j]>1):
                rhog_eq, Tf = isoS_expansion_solver(T0[i,j], rho0_g[i,j], Patm)
                U_boom[i,j] = rhog_eq*EOS_intenergy(rhog_eq, Tf) -U0[i,j]

    #print("U_boom", U_boom)
    U_boom = np.nan_to_num(U_boom )
    return U_boom*1000  #EOS calcuates in kJ


def boom_irr_monophase(rho0g, rho0l, S0_g, S0_l, T0):
    v0   = 1/(S0_l*rho0l + S0_g*rho0g)

    U0 = S0_g*EOS_intenergy(rho0g, T0) + S0_l*EOS_intenergy(rho0l, T0)
    Patm = 101.325
    epsi = 1.0e-4
    deltaxnorm = 10000000.0
    rho = 1.0
    T   = 500.0 #initial guess
    f = np.zeros(2)
    J = np.zeros((2,2))
    Jinv = np.zeros((2,2))
    while(deltaxnorm > 1.0e-4):
      f[0] = boom_f1(rho, T, U0, v0)
      f[1] = EOS_pressure(rho,T)-Patm
      J[0,0] = 0.5*(boom_f1(rho+epsi,T,U0,v0)-boom_f1(rho-epsi,T,U0,v0))/epsi
      J[1,0] = 0.5*(boom_f1(rho,T+epsi,U0,v0)-boom_f1(rho,T-epsi,U0,v0))/epsi
      J[0,1] = 0.5*(EOS_pressure(rho+epsi,T)-EOS_pressure(rho-epsi,T))/epsi
      J[1,1] = 0.5*(EOS_pressure(rho,T+epsi)-EOS_pressure(rho,T-epsi))/epsi
      detJ = J[0,0]*J[1,1]-J[1,0]*J[0,1]
      Jinv[0,0] =  J[1,1]/detJ
      Jinv[1,0] = -J[0,1]/detJ
      Jinv[0,1] = -J[1,0]/detJ
      Jinv[1,1] =  J[0,0]/detJ
      dx    = np.matmul(Jinv,f)
      rho = rho+0.01*dx[0]
      T = T+0.01*dx[1]
      deltaxnorm = np.linalg.norm(dx)
    released_energy = EOS_intenergy(rho,T)-U0
    return released_energy


def boom_f1(rho, T, U0, v0):
    patm = 101.325 
    v = 1.0/rho
    f1 = EOS_intenergy(rho,T)-U0+patm*(v-v0)
    return f1

def boom_irr_biphase(rho0g, rho0l, S0_g, S0_l, T0):
   Patm = 101.325 
   
   v0   = 1/(S0_l*rho0l + S0_g*rho0g)
   U0 = S0_g*EOS_intenergy(rho0g, T0) + S0_l*EOS_intenergy(rho0l, T0)

   if (S0_g<0 or S0_g > 1):
       P0 = EOS_pressure(rho0l,T0)
   else:
       P0 = EOS_pressure(rho0g,T0)

   U0 = EOS_intenergy(rho0,T0)
   rhog_final = 0.5973
   rhol_final = 958.378
   T_final    = 373.12

   vg_final = 1.0/rhog_final
   vl_final = 1.0/rhol_final
   Ug_final = EOS_intenergy(rhog_final, Tfinal)
   Ul_final = EOS_intenergy(rhol_final, Tfinal)
   #f = liquid mass percentage
   #eqn from thiery & mercury 2009
   f = (Patm*( v0-vg_final)-Ug_final+U0)/(Ul_final-Ug_final+Patm*(vl_final-vg_final))
   if ((f<0.0) or (f>1.0)):
     released_energy = np.nan
   else: 
     released_energy = f*Ul_final +(1.0-f)*Ug_final-U0
   return released_energy


def isoS_expansion_solver(T0, rho0, Pout):
    """
    find rho_f, T_f s.t.
    S(rho_f, T_f ) == S(rho0, T0)
    P(rho_f, T_f) == Patm
    """
    S0 = EOS_entropy(rho0,T0)           #
    if (rho0< 100.0):                #       
        rhoout = 1                      #
        Touit = 373                     # initial guesses
    else:                               #
        rhoout = 1000                   #   
        Tout = 373                      #
    x = np.zeros(2)
    x[0] = rho0
    x[1] = T0 
    f = np.zeros(2)
    J = np.zeros((2,2))
    h = 1.0e-4 

    xnorm = 10000000
    oldxnorm = xnorm
    while (xnorm > 1):
        f[0] =  EOS_entropy(x[0], x[1])- S0 
        f[1] = EOS_pressure(x[0], x[1]) - Pout
        J[0,0] = (EOS_entropy(x[0]+h, x[1])-EOS_entropy(x[0]-h, x[1]))/(2*h)
        J[0,1] = (EOS_entropy(x[0], x[1]+h)-EOS_entropy(x[0], x[1]-h))/(2*h)
        J[1,0] = (EOS_pressure(x[0]+h, x[1])-EOS_pressure(x[0]-h, x[1]))/(2*h)
        J[1,1] = (EOS_pressure(x[0], x[1]+h)-EOS_pressure(x[0], x[1]-h))/(2*h)

        Jinv = np.linalg.inv(J)
        dx = -Jinv@f 
        oldxnorm = xnorm
        xnorm = np.linalg.norm(dx)
        #dprint(oldxnorm, xnorm)
        if (xnorm > 10*oldxnorm):
            x = x +0.01*dx
        elif (xnorm > 100*oldxnorm):
            x = x +0.001*dx
        else:
            x = x+0.1*dx
    return x[0], x[1]



   
def validate():  
   T = 647
   rho = 358
   delta = rho/rho_c
   tau = Tc/T
   assert (abs(get_phi_residual(delta, tau)-(-0.121202657e+1))<1e-6)
   assert (abs(get_dPhiR_tau(delta, tau)-(-0.321722501e+1))<1e-6)
   assert (abs(get_dPhiR_delta(delta, tau)-(-0.714012024))<1e-6)
   assert (abs(get_phi_ideal(delta, tau)-(-0.156319605e+1))<1e-6)
   assert (abs(get_dPhi0_tau(delta, tau)-(0.980343918e+1))<1e-6)
   assert (abs(EOS_pressure(996.513, 300)-3.537< 1e-6))
   assert (abs(EOS_pressure(0.02559, 300)-3.537< 1e-6))
   assert (abs(EOS_pressure(0.42559, 1273)-250< 1e-6))

validate()
