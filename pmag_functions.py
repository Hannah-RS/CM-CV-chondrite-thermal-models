"""Functions to calculate a quantity proportional to pTRM of a sample"""
import numpy as np
from scipy.interpolate import interp1d

#interpolate function for TRM acquisition from paleomagnetic data
#magnetite
magdata = np.loadtxt('../CV_chondrites/ptrm-mag.csv',skiprows=1, delimiter=',')
ptrm_func_mag = interp1d(magdata[:, 0], magdata[:, 1])
#pyrrhotite
pyrrhdata = np.loadtxt('../CV_chondrites/ptrm-pyrrh.csv',skiprows=1, delimiter=',')
ptrm_func_pyrrh = interp1d(pyrrhdata[:, 0], pyrrhdata[:, 1])
#taenite
taedata = np.loadtxt('../CV_chondrites/ptrm-tae.csv',skiprows=1, delimiter=',')
ptrm_func_tae = interp1d(taedata[:, 0], taedata[:, 1])

def nrm(Tblock,Tmub,Tneb,carrier,carrier_min):
    """Calculate a proxy for the NRM of a carrier within a sample.
        Currently assumes a linear relationship.
        Parameters
        ----------
        Tblock : float
            Blocking temperature of the carrier (K)
        Tmub : float
            Maximum unblocking temperature of the sample (K)
        Tneb : float
            Temperature of the sample at the time of nebula field dissipation (K)
        carrier : str
            Type of carrier (magnetite, pyrrhotite, taenite)
        carrier_min : float
            Minimum temperature of the experimental data
        Returns
        -------
        nrm : float
            Proxy for the NRM of the sample"""
    if (Tmub >= Tblock) & (Tneb < Tblock): #Fully magnetized
        ptrm_mub = 0 
    elif (Tmub >= Tblock) & (Tneb >= Tblock): #Fully unmagnetized
        ptrm_mub = 0 
    elif (Tmub < Tblock): #Partially magnetized
        if carrier == 'magn':
            ptrm_mub = ptrm_func_mag(Tmub)
        elif carrier == 'pyrrh':
            ptrm_mub = ptrm_func_pyrrh(Tmub)
        elif carrier == 'tae':
            ptrm_mub = ptrm_func_tae(Tmub)
    if (Tneb > carrier_min) & (Tneb < Tblock): #Partially unmagnetized, 
        #273 cut off because of interpolation limits
        if carrier == 'magn':
            ptrm_neb = ptrm_func_mag(Tneb)
        elif carrier == 'pyrrh':
            ptrm_neb = ptrm_func_pyrrh(Tneb)
        elif carrier == 'tae':
            ptrm_neb = ptrm_func_tae(Tneb)
    else:
        ptrm_neb = 1
    nrm = ptrm_neb - ptrm_mub
    return nrm

def multi_carry(nrma,f):
    """Calculate a proxy for the NRM of a sample with two carriers.
        Parameters
        ----------
        nrma : np.ndarray
            array of values of proxy for the NRM of each carrier
        f : np.ndarray
            array of values of the fraction of each carrier 
        Returns
        -------
        nrmt : float
            Proxy for the total NRM of the sample"""
    nrmt = np.sum(nrma*f)
    return nrmt