"""Plot nrm as a function of accretion time and size."""


import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma


cmap = 'viridis'
norm = mpl.colors.Normalize(vmin=0, vmax=1)

#load saved results
sts = np.load('../Results/sts.npz')['data']
ss = np.load('../Results/ss.npz')['data']
nrmtot = np.load('../Results/nrmtot.npz')['data']

#Calculate the max values but mask the nans for plotting.
fig, ax = plt.subplots(1,3,sharey='row',figsize=(8,3),tight_layout=True)
ax[0].pcolormesh(sts,ss/1000,np.transpose(np.min(nrmtot,axis=2)),cmap=cmap, norm=norm)
ax[1].pcolormesh(sts,ss/1000,np.transpose(np.average(nrmtot,axis=2)),cmap=cmap, norm=norm)
ax[2].pcolormesh(sts,ss/1000,np.transpose(np.nanmax(nrmtot,axis=2)),cmap=cmap, norm=norm)
cax = fig.add_axes([1.01, 0.2, 0.02, 0.5])
cb=plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm),cax=cax,orientation='vertical')

ax[0].set(ylabel='Planetesimal radius (km)', title='minimum pTRM',
          xlabel='Accretion time (Myr)',xlim=[1.6,5])
ax[1].set(title='average pTRM',xlabel='Accretion time (Myr)',
          xlim=[1.6,5])
ax[2].set(title='maximum pTRM',xlabel='Accretion time (Myr)',
          xlim=[1.6,5])
fig.suptitle('all CV chondrites')
plt.savefig('../Plots/nrm.png', dpi=300, bbox_inches='tight')