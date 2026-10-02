"""Functions to calculate a quantity proportional to pTRM of a sample"""
import numpy as np
from scipy.interpolate import interp1d
from constants import fTblock_min

#interpolate function for TRM acquisition from paleomagnetic data
pmag_data =np.loadtxt('../CV_chondrites/Kaba_TRM.csv',skiprows=4,delimiter=',')
#normalised pTRM as a function of T/Tblock in K
ptrm_func = interp1d(pmag_data[:,2],pmag_data[:,7])

def nrm(Tblock,Tmub,Tneb):
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
        Returns
        -------
        nrm : float
            Proxy for the NRM of the sample"""
    #pTRM at the maximum unblocking temperature
    if (Tmub >= Tblock): #Could be fully magnetized
        ptrm_mub = 0 
    elif (Tmub < Tblock): #Partially magnetized
        ptrm_mub = ptrm_func(Tmub/Tblock)

    #pTRM at the time of nebula field dissipation
    if (Tneb >= Tblock): #Fully unmagnetized
        ptrm_neb = 0
    elif (Tneb < Tblock) & (Tneb/Tblock > fTblock_min): #Partially unmagnetized
        ptrm_neb = ptrm_func(Tneb/Tblock)
    elif (Tneb < Tblock) & (Tneb/Tblock <= fTblock_min): #At data limit, assume = 1
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