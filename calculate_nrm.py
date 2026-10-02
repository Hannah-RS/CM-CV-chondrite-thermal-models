""" Calculate NRM of each meteorite for a given accretion time and planetesimal size. Save the output."""
import numpy as np
import pandas as pd
from pmag_functions import nrm, multi_carry
from constants import tbpyrrh, tbmag, tbtae, ns, nst#, CVtype

CVtype = 'Red'
#load savedtemperature at the time of nebula field dissipation
Tneb = np.load(f'Results/{CVtype}/Tneb.npz')['data']

#load Clara data
cv_dat = pd.read_csv('../CV_chondrites/Table-for-Hannah-0926.csv')
cv_dat['max T (K)'] = cv_dat['Max T (deg C)'] + 273
if CVtype != 'all':
    cv_dat = cv_dat[cv_dat['Group'] == CVtype] #filter for subset of CV chondrites
nm = len(cv_dat['Meteorite']) #number of meteorites

nrmtot = np.zeros((nst,ns,nm)) #NRM of each meteorite for each accretion time and planetesimal size
Tblocks = np.array([tbpyrrh, tbmag, tbtae]) #blocking temperatures of each carrier
carriers = ['pyrrh', 'magn', 'tae'] #names of each carrier
for i in range(nst):
    for j in range(ns):
        for k, Tmub in enumerate(cv_dat['max T (K)']):
            nrma = np.zeros(3) #array to store NRM of each carrier
            f = cv_dat.iloc[k,5:8].to_numpy() #fraction of each carrier
            if np.isnan(Tneb[i,j,k]): #if Tneb is nan, set NRM to nan
                nrmtot[i,j,k] = np.nan
            else:
                for Tblock, carrier in zip(Tblocks,carriers):
                    nrma = nrm(Tblock, Tmub, Tneb[i,j,k])
                nrmtot[i,j,k] = multi_carry(nrma, f)

np.savez(f'Results/{CVtype}/nrmtot.npz', data=nrmtot)