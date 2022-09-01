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
stress2=np.loadtxt('OUT\pi_plane2.csv', delimiter=',', dtype=np.float64)
dev_h=np.loadtxt('OUT\DEV_H.csv', delimiter=',', dtype=np.float64)
#dev_s=np.loadtxt('OUT\DEV_S.csv', delimiter=',', dtype=np.float64)
dev_s=np.zeros([3])
dev_s[1]=1.0
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
#########################################################################################
# PLOT PI-PLANE & WRITE OUTPUT DATA
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
plt.figure(1,figsize=(8,8))
#plt.title('$\pi$-Plane', size=20)
plt.plot([0,pp_axis[0,0]],[0,pp_axis[0,1]], linewidth=3, color='k')
plt.annotate("$s_{xx}$", xy=(pp_axis[0,0],pp_axis[0,1]), size=30, xytext=(pp_axis[0,0]-0.0,pp_axis[0,1]-0.15))
plt.plot([0,pp_axis[1,0]],[0,pp_axis[1,1]], linewidth=3, color='k')
plt.annotate("$s_{yy}$", xy=(pp_axis[1,0],pp_axis[1,1]), size=30, xytext=(pp_axis[1,0]+0.05,pp_axis[1,1]-0.05))
plt.plot([0,pp_axis[2,0]],[0,pp_axis[2,1]], linewidth=3, color='k')
plt.annotate("$s_{zz}$", xy=(pp_axis[2,0],pp_axis[2,1]), size=30, xytext=(pp_axis[2,0]+0.05,pp_axis[2,1]-0.05))
ang=360;
#dev_s=dev_s*0.9;
dev_s=dev_s*0.575;
#dev_s=dev_s*0.25;
plt.arrow(0,0,dev_s[0],dev_s[1], width=0.015, color='tab:purple', zorder=5)
dev_h=dev_h*0.5;
plt.arrow(0,0,dev_h[0],dev_h[1], width=0.015, color='tab:red', zorder=5)
r2=r*1.1
plt.axis([-r2,r2,-r2,r2])
plt.axis('off')

# Original yield surface
plt.plot(stress[:,2], stress[:,3], c='gray', linestyle='--', linewidth=2.5, label='Isotropic hardening')
# HAH20h
plt.plot(stress[:,0], stress[:,1], c='tab:blue', linewidth=2.5, label='HAH'+'$_{20h}$')
# HAH20e
#plt.plot(stress2[:,0], stress2[:,1], linewidth=2.5, label='HAH'+'$_{20\epsilon}$')
#plt.legend()
plt.show()


