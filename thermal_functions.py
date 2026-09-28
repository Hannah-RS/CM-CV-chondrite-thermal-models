"""Functions for calculating the thermal history of a chondritic body taken from the Jupyter Notebook in the Sanderson_et_al_2026 directory"""

import numpy as np
from constants import LAs, H0s, SSAGE, LN2

def calcOriginAbund(abund,el):
    '''
    Given a present-day abundance of a radioisotope, calculates the abundance at the beginning of the solar system.
    '''
    if el not in LAs:
        print("Element "+el+" not in radiosotope list.")
    LAel=LAs[el]
    HLs=SSAGE/LAel
    return abund*(2**HLs)


def calcHeatProd(el,conc, cp, dt):
    '''
    Calculates the heat produced by the decay of a radioisotope
        el: radioisotope (string, one of the isotopes in the above dictionaries)
        conc: Current concentration of the isotope (fraction)
        cp: heat capacity (J/kg/K)
        dt: time step (seconds)
    '''
    return (H0s[el]*conc/cp*dt)

def calcDecayFrac(el,C0,time):
    '''
        el: radioisotope (string, one of the isotopes in the above dictionaries)
        C0: concentration (fraction)
        time: (seconds)
    '''
    return C0*np.exp(-LN2*time/LAs[el])


def cell_decay(dtime,sc):
    """
       Updates the abundances of Al26 after decaying for a set interval of time, dtime.
    """
    return calcDecayFrac('Al26',sc,dtime)


def cell_heat(sc,dtime,Cp):
    """
       Calculates the heat produced from decaying aluminium-26
    """
    DT = 0   # change in temperature
    i='Al26'
    DT = DT + calcHeatProd(i, sc, Cp, dtime)

    return DT

def cell_heat_K(sc,dtime,Cp):
    """
       Calculates the heat produced from decaying potassium-40.
    """
    DT = 0   # change in temperature
    i='K40'
    DT = DT + calcHeatProd(i, sc, Cp, dtime)

    return DT


def calcStep(Cp, K, rho, DR, timeC):
    """
       Calculates minimum time step.
    """
    alpha = K/(Cp*rho)
    DT = np.min(DR**2/alpha/3)
    return DT

def thermalCondStep(T,Cp,K,rho,Rads,dt):
    """
       Thermal coduction finite difference approximation.
    """
    nc = len(T)

    C = Rads
    DR = np.diff(C,append=0)
    DR[-1]=DR[-2]

    DT = dt

    D=np.zeros((nc,3))

    for i in range(0,nc):
        for j in range(0,3):
            if i==0:
                D[0,j]=K[0]/(rho[i]*Cp[0])
            elif i==nc-1:
                D[nc-1,j]=K[nc-1]/(rho[i]*Cp[nc-1])
            else:
                D[i,j]=(K[i]+K[i+(j-1)])/(rho[i]*Cp[i])/2.0

    dn = D[:,1]
    dp = D[:,2]
    dm = D[:,0]

    diagU = np.zeros(nc)
    diagC = np.zeros(nc)
    diagL = np.zeros(nc)

    for i in range(1,nc-1):
        diagL[i] = DT/(2.0*DR[i]*DR[i+1])*(dn[i]/(i)-dm[i])
        diagC[i] = 0.5+DT/(2.0*DR[i]*DR[i+1])*(dp[i]+dm[i])
        diagU[i] = -DT/(2.0*DR[i]*DR[i+1])*(dn[i]/(i)+dp[i])

    diagC[0] = 0.5 + DT/(DR[0]**2)*dn[0]
    diagU[0] = -DT/(DR[0]**2)*dn[0]

    # Calc source
    Hadd = np.zeros(nc)
    for i in range(0,nc-1):
        Hadd[i] = 0.5*T[i]
    Hadd[nc-2]=Hadd[nc-2]+DT/(2*DR[nc-2]**2)*(dn[nc-2]/(nc-2)+dp[nc-2])*T[nc-1]

    # Solve
    U = np.zeros(nc)
    GAM = np.zeros(nc)

    BET=diagC[0]
    U[0] = Hadd[0]/BET
    U[-1]=T[-1]

    for i in range(1,nc-1):
        GAM[i]=diagU[i-1]/BET
        BET=diagC[i]-diagL[i]*GAM[i]
        U[i] = (Hadd[i]-diagL[i]*U[i-1])/BET
    for i in range(nc-2,-1,-1):
        U[i]=U[i]-GAM[i+1]*U[i+1]

    return U


def calc_Cp(temp,WR):
    """
       Calculates Cp of rock and ice mix.
    """
    WR = WR*(temp<273)
    
    Cp_ice = 185+7.037*(temp)
    Cp_rock= 920
    
    Cp=Cp_ice*WR/(WR+1)+Cp_rock/(WR+1)
    
    return Cp
    
def calc_k(temp,WR):
    """
       Calculates thermal conductivity of rock and ice mix.
    """
    WR = WR*(temp<273)
    
    k_ice = 0.4685+488.12/(temp)
    k_rock= 4.2
    
    k=k_ice*WR/(WR+1)+k_rock/(WR+1)
    
    return k
    
def calc_rho(temp,WR):
    """
       Calculates density of rock and ice mix.
    """
    r_ice = 917
    r_rock= 3000
    
    r=r_ice*WR/(WR+1)+r_rock/(WR+1)
    
    return r+0*temp

def add_latent_heat(temp, new_temp, cps, WR):
    """
       adds latent heat.
    """
    ice_melt = 3.34e5  #J/kg
    rock_hyd = 2.875e5 #J/kg
    
    M_ice=WR/(WR+1)
    M_rock=1/(WR+1)
    
    Mf_react = 240/36
    M_rock_hyd_max=Mf_react*M_ice
    
    if M_rock_hyd_max<M_rock:
        M_rock=M_rock_hyd_max
        
    defecit=np.zeros(len(cps))
        
    for i in range(0,len(temp)):
        if temp[i]<273 and new_temp[i]>273:
            dT=(M_rock*rock_hyd/cps[i]-M_ice*ice_melt/cps[i])
            if dT<0:
                defecit[i]=dT
                dT=0
            new_temp[i]=new_temp[i]+dT
    
    return new_temp, defecit

def timeStep(temp, WR, nr, Rs, dt, Al26, defecit):
    """
       calculates thermal properties, adds decay heat, executes thermal conduction step, adds latent heat.
    """
    cps=calc_Cp(temp,WR)
    ks=calc_k(temp,WR)
    rhos=calc_rho(temp,WR)
    
    for i in range(0,nr):
        if cps[i]<500 or cps[i]>5000:
            print("Cp "+str(cps[i]))
        if ks[i]<0.5 or ks[i]>5:
            print("k "+str(ks[i]))
        if rhos[i]<500 or rhos[i]>5000:
            print("rho "+str(rhos[i]))
        
    for i in range(0,nr):
        # do not add heat to cell on surface
        if not i==nr-1:
            temp[i]+=cell_heat(Al26,dt,cps[i])
            temp[i]+=cell_heat_K(700e-9,dt,cps[i])
    
    new_temp=thermalCondStep(temp,cps,ks,rhos,Rs,dt)
    
    dt=new_temp-temp
    for i in range(0,nr):
        if defecit[i]>0 and dt[i]>0:
            if dt[i]>defecit[i]:
                new_temp[i]=new_temp[i]-defecit[i]
                defecit[i]=0
            else:
                new_temp[i]=temp[i]
                defecit[i]=defecit[i]-dt[i]
    
    new_temp, new_defecit=add_latent_heat(temp, new_temp, cps, WR)
    
    defecit=defecit+new_defecit
    
    return new_temp, defecit


def runModel(start, nt, dt, Rs, nr, rad, sT, Al26, WR):
    """
       Executes model
    """
    
    #cps = np.ones((nt,nr))*Cp
    #rhos = np.ones((nt,nr))*rho
    #ks = np.ones((nt,nr))*k
    
    Al26_arr = np.zeros(nt)
    
    dr = rad/nr
    
    temp = np.ones((nt,nr))*sT
    
    Al26=cell_decay(start,Al26)
    
    Al26_arr[0]=Al26
    
    defecit=np.zeros(nr)
    for i in range(0,nt-1):
        temp[i+1,:],deficit=timeStep(temp[i,:], WR, nr, Rs, dt, Al26, defecit)
        #temp[i+1,:]=timeStep(temp[i,:], ks[i,:], rhos[i,:], cps[i,:], nr, Rs, dt, Al26)
        Al26=cell_decay(dt,Al26)
        Al26_arr[i]=Al26
        
        #ks[i+1,:]=calc_k(temp[i+1,:])
        #rhos[i+1,:]=calc_Cp(temp[i+1,:])
        #cps[i+1,:]=calc_rho(temp[i+1,:])
        
    return temp, Al26_arr
        
def calc_shell_mass(r,dr,rho):
    """
       Calculates mass of spherical shell.
    """
    V= 4/3 * np.pi * (r**3 - (r-dr)**3) 
    M= V*rho
    return M

def mass_temp_hist(R,dr,nr,rho,T):
    """
       Generates a histogram of mass by temperature for thermal model.
    """
    M = np.zeros(nr)
    for i in range(0,nr):
        M[i] = calc_shell_mass(R[i],dr,rho)
    
    #fig, ax = plt.subplots()
    #ax.hist(T, nbins, weights=M)
    
    nbins=np.linspace(100,2000,30)
    return np.histogram(T, nbins, weights = M)

