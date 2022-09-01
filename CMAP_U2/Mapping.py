import math
import numpy as np
import numpy.linalg as lin
import pandas as pd
from matplotlib import cm
from matplotlib import ticker
#from UMAT import *
import matplotlib.pyplot as plt
import csv
import sys

############################################################################
#                                                                          #
#                               MAIN PROGRAM                               #
#                                                                          #
############################################################################
global hard_par, yld_par, flow_par, sua_par
global qm, pm, k, k1, k2, k3, k4, k5
global c, kc, l, kl
global am, a1, a2, a3 ,a4, a5, a6, a7, a8
global p1, p2, p3, p4, p5, p6, p7, p8, p9
global ndim1, ndim2, ndim3, ndim4, ndim5, ndim6, ndim7
# SET PLOTTING PARAMETERS #
print_par=1       # Print the result on csv file. (0: off / 1: on)
r=1.2             # RADIUS OF AXIS
# OPEN DATA FILES #
stress=np.loadtxt('OUT\pi_plane.csv', delimiter=',', dtype=np.float64)
cmap_data = np.loadtxt('OUT\CMAP.csv', delimiter=',', dtype=np.float64)
debug = np.loadtxt('DEBUG\DEBUG.csv', delimiter=',', dtype=np.float64)
indx=0
ndiv=0
nconv=0
#ang1=721
#ang2=361
ang1=max(debug[:,0])-min(debug[:,0])+1
ang2=max(debug[:,1])-min(debug[:,1])+1
#sys.exit()
ndata=ang1*ang2
ang1=int(ang1)
ang2=int(ang2)
cmap_output=np.zeros([ang1,ang2,3])
cmap_fail=np.zeros([ang1,ang2,3])
for i in range(ang1):
      for j in range(ang2):
            if cmap_data[indx,2] == 200:
                  cmap_fail[i,j,0]= cmap_data[indx,0]
                  cmap_fail[i,j,1]= cmap_data[indx,1]
                  cmap_fail[i,j,2]= cmap_data[indx,2]
                  ndiv=ndiv+1
            else:
                  cmap_output[i,j,0]=cmap_data[indx,0]
                  cmap_output[i,j,1]=cmap_data[indx,1]
                  cmap_output[i,j,2]=cmap_data[indx,2]
                  nconv=nconv+1
            indx=indx+1
#########################################################################################
# PLOT PI-PLANE & WRITE OUTPUT DATA
# OUTLINE OF PI-PLANE
pp_axis=np.zeros([3,2])
th=0*np.arccos(-1)/180
pp_axis[0,0]=r/np.sqrt(1+np.tan(th)**2)
pp_axis[0,1]=np.tan(th)*pp_axis[0,0]
th=120*np.arccos(-1)/180
pp_axis[1,0]=-r/np.sqrt(1+np.tan(th)**2)
pp_axis[1,1]=np.tan(th)*pp_axis[1,0]
th=240*np.arccos(-1)/180
pp_axis[2,0]=-r/np.sqrt(1+np.tan(th)**2)
pp_axis[2,1]=np.tan(th)*pp_axis[2,0]
plt.figure(1,figsize=(10,8))
#plt.figure(1,figsize=(8,8))
#plt.title('$\pi$-Plane', size=20)
#plt.arrow(0,0,pp_axis[0,0],pp_axis[0,1], width=0.01, color='black')
plt.plot([0,pp_axis[0,0]],[0,pp_axis[0,1]], linewidth=3, color='k')
plt.annotate("$s_{xx}$", xy=(pp_axis[0,0],pp_axis[0,1]), size=30, xytext=(pp_axis[0,0]-0.1,pp_axis[0,1]-0.15))
#plt.arrow(0,0,pp_axis[1,0],pp_axis[1,1], width=0.01, color='black')
plt.plot([0,pp_axis[1,0]],[0,pp_axis[1,1]], linewidth=3, color='k')
plt.annotate("$s_{yy}$", xy=(pp_axis[1,0],pp_axis[1,1]), size=30, xytext=(pp_axis[1,0]+0.05,pp_axis[1,1]-0.05))
#plt.arrow(0,0,pp_axis[2,0],pp_axis[2,1], width=0.01, color='black')
plt.plot([0,pp_axis[2,0]],[0,pp_axis[2,1]], linewidth=3, color='k')
plt.annotate("$s_{zz}$", xy=(pp_axis[2,0],pp_axis[2,1]), size=30, xytext=(pp_axis[2,0]+0.05,pp_axis[2,1]-0.05))
r2=r*1.02
plt.axis([-r2,r2,-r2,r2])
plt.axis('off')
#max_iter=np.amax(cmap_output[:,:,2])
nfrac1=ndiv/ndata*100.0
nfrac2=nconv/ndata*100.0
print(ndata, nconv)
mean_iter=np.mean(cmap_output[:,:,2])
#mean_iter=nfrac2
props=dict(boxstyle='square', facecolor='wheat', alpha=0.5)
plt.text(0.2, 1.0, 'Mean iter: '+str(int(mean_iter)),size=25, bbox=props)
#props2=dict(boxstyle='square', facecolor='white', alpha=0.5)
#lt.text(-1.5, 1.0, r'(a)',size=30, bbox=props2)
plt.plot(stress[:,0], stress[:,1], 'k', linewidth=2.5)
plt.plot(stress[:,2], stress[:,3], 'k--', linewidth=2.5)

max_iter=5 # CTM
step=6
#max_iter=20 # CPPM
#max_iter=40 # LSEBM

plt.scatter(cmap_output[:,:,0], cmap_output[:,:,1], c=cmap_output[:,:,2], cmap='Blues', s=2, vmin=0, vmax=max_iter)

cb=plt.colorbar()
cb.ax.set_yticklabels(labels=np.linspace(0,max_iter,step), size=25)
cb.locator = ticker.MaxNLocator(nbins=step)
cb.update_ticks()

plt.scatter(cmap_fail[:,:,0], cmap_fail[:,:,1], c='lightcoral', s=2)

plt.show()

