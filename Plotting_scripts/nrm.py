"""Plot nrm as a function of accretion time and size."""


import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma

CVtypes = ['all', 'OxA', 'OxB', 'Red']

cmap = 'viridis'
norm = mpl.colors.Normalize(vmin=0, vmax=1)



#Calculate the max values but mask the nans for plotting.
fig, ax = plt.subplots(4,3,sharey='row',sharex='col',figsize=(8,9),tight_layout=True)
for i, CVtype in enumerate(CVtypes):
    #load saved results
    sts = np.load(f'../Results/{CVtype}/sts.npz')['data']
    ss = np.load(f'../Results/{CVtype}/ss.npz')['data']
    nrmtot = np.load(f'../Results/{CVtype}/nrmtot.npz')['data']
    #plot
    ax[i,0].pcolormesh(sts,ss/1000,np.transpose(np.min(nrmtot,axis=2)),cmap=cmap, norm=norm)
    ax[i,1].pcolormesh(sts,ss/1000,np.transpose(np.average(nrmtot,axis=2)),cmap=cmap, norm=norm)
    ax[i,2].pcolormesh(sts,ss/1000,np.transpose(np.nanmax(nrmtot,axis=2)),cmap=cmap, norm=norm)
    ax[i,0].set(ylabel='Planetesimal radius (km)', title=f'minimum {CVtype}')
    ax[i,1].set(title=f'average {CVtype}')
    ax[i,2].set(title=f'maximum {CVtype}')
for i in range(3):
    ax[3,i].set(xlabel='Accretion time (Myr)',xlim=[1.6,2.8])
cax = fig.add_axes([1.01, 0.2, 0.02, 0.5])
cb=plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm),cax=cax,orientation='vertical',label='pTRM')

#plt.savefig(f'../Plots/nrm_{CVtype}.png', dpi=300, bbox_inches='tight')