"""Constants required for calculating the thermal history of a chondritic body taken from the Jupyter Notebook in the Sanderson_et_al_2026 directory"""

import numpy as np

### CONSTANTS
YR = 3.1536E7             # converting a year to seconds
LN2 = np.log(2)
SSAGE=4.571*10**9*YR      # age of solar system


dH_ice=3.34e5             # latent heat of ice melting (J/kg)
dH_hyd_mol=69e3           # heat of hydration (J/mol)
dH_hyd=dH_hyd_mol/240*1e3 # Coverting to J/kg rock altered

Mf_react = 36/240              
melt_vs_react = dH_ice/dH_hyd 
WR = Mf_react*melt_vs_react   

tstart = 1.6
tstop = 5
nst=18 # Number of discrete accretion times.
ns=61 #Number of discrete radii.
rstart = 5000
rstop = 150000

# Radioactive isotopes:
# Heat Production (W/kg)
H0U238  = 94.65E-6
H0U235  = 568.7E-6
H0TH    = 26.38E-6
H0K     = 29.17E-6
H0FE    = 6.7E-2
H0AL    = 3.57E-1
H0MN    = 5.E-3
H0CA    = 1.53E-1

# Half-Lives
LAU238  = 4.47E9*YR
LAU235  = 7.0381E8*YR
LATH    = 14.01E9*YR
LAK     = 1.277E9*YR
LAFE    = 1.5E6*YR
LAAL    = 7.16E5*YR
LAMN    = 3.74E6*YR
LACA    = 0.103E5*YR

# Put in dictionary
LAs = dict({'U238': LAU238, 
              'U235':LAU235,
              'Th232':LATH,
              'K40':LAK,
              'Fe60':LAFE,
              'Al26':LAAL,
              'Mn53':LAMN,
              'Ca41':LACA})
H0s = dict({'U238': H0U238, 
              'U235':H0U235,
              'Th232':H0TH,
              'K40':H0K,
              'Fe60':H0FE,
              'Al26':H0AL,
              'Mn53':H0MN,
              'Ca41':H0CA})


#Blocking temperatures from Tauxe Essentials of Paleomagnetism and Nichols 2021
tbpyrrh = 325 + 273 #pyrrhotite blocking temperature (K)
tbmag = 580 + 273 #magnetite blocking temperature (K)
tbtae = 360 + 273 #taenite blocking temperature (K)
tneb = 4.5*1e6*YR #time of dissipation of the nebula field (s)