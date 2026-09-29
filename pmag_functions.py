"""Functions to calculate a quantity proportional to pTRM of a sample"""
import numpy as np

@np.vectorize
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
    if (Tmub >= Tblock) & (Tneb < Tblock): #Fully magnetized
        nrm = 1 
    elif (Tmub >= Tblock) & (Tneb >= Tblock): #Fully unmagnetized
        nrm = 0 
    elif (Tmub < Tblock): #Partially magnetized
        nrm = (Tmub - Tneb)/Tblock     #This relationship will probably get changed at some point.
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