"""Contour plot of time cooled through pyrrhotite blocking temperature for locations that are consistent with CV chondrite peak metamorphic temperatures"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


cmap_bin = matplotlib.colors.ListedColormap(['green', 'grey'])

#load saved results
sts = np.load('../Results/sts.npz')['data']
ss = np.load('../Results/ss.npz')['data']
tpyrrh = np.load('../Results/tpyrrh.npz')['data']

plt.figure(figsize=(8,6))
plt.contourf(sts,ss/1000,np.transpose(tpyrrh),cmap=cmap_bin, levels=[-0.5,4.5,100])
cb=plt.colorbar()
cb.set_ticks([3,50], labels=['Nebula field active','Nebula field dissipated'])
cb.ax.set_ylim(0,100)

CS=plt.contour(sts,ss/1000,np.transpose(tpyrrh),np.linspace(0,10,5),colors='w',linewidths=0.5)
plt.clabel(CS, inline=1, fontsize=12)
plt.ylabel('Planetesimal radius (km)')
plt.xlabel('Accretion time (Myr)')

plt.xlim([1.6,5])
plt.title('Time cooled through pyrrhotite blocking temperature (593 K)')
plt.savefig('../Plots/CV_block.png', dpi=300, bbox_inches='tight')