"""Plot maximum depth of C V chondrites"""

import matplotlib.pyplot as plt
import numpy as np

#load saved results
sts = np.load('../Results/sts.npz')['data']
ss = np.load('../Results/ss.npz')['data']
max_depth = np.load('../Results/max_depth.npz')['data']
min_depth = np.load('../Results/min_depth.npz')['data']

fig, ax = plt.subplots(1,2,sharey='row',tight_layout=True,figsize=(12,6))
vmin=1
vmax = 50
cmap_depth= plt.cm.viridis
cmap_depth.set_under('lightgray')
im0 = ax[0].pcolormesh(sts, ss / 1000, np.transpose(max_depth) / 1e3, cmap=cmap_depth, vmin=vmin, vmax=vmax)
im1 = ax[1].pcolormesh(sts, ss / 1000, np.transpose(min_depth) / 1e3, cmap=cmap_depth, vmin=vmin, vmax=vmax)
cax = fig.add_axes([1, 0.15, 0.02, 0.7])
cbar = fig.colorbar(im0, cax=cax)
cbar.set_label('Depth (km)')
#plt.colorbar(cmap=plt.cm.viridis, ax=ax[0], label='Depth (km)')


ax[0].set(ylabel='Planetesimal radius (km)',xlabel='Accretion time (Myr)',title='Maximum depth of CV chondrites (km)',xlim=[1.6,5])
ax[1].set(xlabel='Accretion time (Myr)',title='Minimum depth of CV chondrites (km)',xlim=[1.6,5])
plt.savefig('../Plots/CV_chondrite_depth.png', dpi=300,bbox_inches='tight')
