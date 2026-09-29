"""Run thermal models for a range of sizes and formation times to compare with CV chondrite peak unblocking temperatures"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
from constants import YR, tneb
from thermal_functions import runModel


#Import data from Clara
cv_dat = pd.read_csv('../CV_chondrites/Summary-TH-data-for-Hannah.csv')
cv_dat['max T (K)'] = cv_dat['Max T (deg C)'] + 273
nm = len(cv_dat['Meteorite']) #number of meteorites

nst=18
sts=np.linspace(1.6,5,nst)

ns=31
ss=np.linspace(5000,150000,ns)

maxT = np.zeros((nst,ns))
CVagree = np.zeros((nst,ns)) #Could the CV chondrite data have come from this body?
depth = np.zeros((nst,ns,nm)) #Depth of each meteorite in the body (m)  
Tneb = np.zeros((nst,ns,nm)) #Temperature at each meteorite position the time of the dissipation of the nebula field (K)

for i in range(0,nst):
    for j in range(0,ns):
        WR=0.7

        nr=100
        rad=ss[j]
        nt=300
        dt=4e4*YR
        sT=180

        Rs=np.linspace(0,rad,nr)
        ts=np.arange(0,nt*dt,dt)

        factor=1
        Al_tot_abund=factor * 1.14/100
        Al_26_27_start=5*10**-5   
        Al26=Al_tot_abund*Al_26_27_start
        
        start=sts[i]*YR*1e6
        temp_map=np.zeros((nt,nr))
        temp_map,Al26_arr=runModel(start, nt, dt, Rs, nr, rad, sT, Al26, WR)
        
        maxT[i,j] = np.max(temp_map)

        #find depth with given max temperature
        Tmax_wdepth = np.max(temp_map,axis=0)
        for Tval in cv_dat['max T (K)']:
            if np.any(Tmax_wdepth > Tval)&(np.all(Tmax_wdepth < 1360)): #Check at least one position exceeds peak temp and no silicate melting
                CVagree[i,j] = 1
            else:
                CVagree[i,j] = 0
                break
        if CVagree[i,j] == 1:
            for k, Tval in enumerate(cv_dat['max T (K)']):
                idr = np.where(Tmax_wdepth >= Tval)[0][-1] #position of location with this peak temp - shallowest position it is hotter
                depth[i,j,k] = rad-Rs[idr]
                idt = np.where(ts>tneb)[0][0] #index of time of nebula field dissipation
                Tneb[i,j,k] = temp_map[idt,idr] #Temperature at the time of nebula field dissipation
Tneb[CVagree==0] = np.nan

#save results to npz files
variables_to_save = {
    'maxT': maxT,
    'CVagree': CVagree,
    'depth': depth,
    'Tneb': Tneb,
    'sts': sts,
    'ss': ss,
}

for name, data in variables_to_save.items():
    np.savez(f'Results/{name}.npz', data=data)

