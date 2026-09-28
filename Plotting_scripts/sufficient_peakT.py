"""Contour plot of accretion times and radii consistent with the CV chondrites"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

#load saved results
CVagree = np.load('../Results/CVagree.npz')['data']
sts = np.load('../Results/sts.npz')['data']
ss = np.load('../Results/ss.npz')['data']
maxT = np.load('../Results/maxT.npz')['data']


cmap_bin = matplotlib.colors.ListedColormap(['grey', 'green'])

plt.figure(figsize=(8,6))
plt.contourf(sts,ss/1000,np.transpose(CVagree),cmap=cmap_bin, levels=[-0.5,0.5,1.5])
cb=plt.colorbar()
cb.set_ticks([0,1], labels=['Inonsistent','Consistent'])
cb.ax.set_ylim(-0.5,1.5)

CS=plt.contour(sts,ss/1000,np.transpose(maxT),np.linspace(0,1400,8),colors='w',linewidths=0.5)
plt.clabel(CS, inline=1, fontsize=12)
plt.ylabel('Planetesimal radius (km)')
plt.xlabel('Accretion time (Myr)')

plt.xlim([1.6,5])
plt.savefig('../Plots/CV_consistent.png', dpi=300, bbox_inches='tight')
