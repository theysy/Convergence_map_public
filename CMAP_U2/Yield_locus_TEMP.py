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
# OPEN DATA FILES #
stress=np.loadtxt('OUT\yld_locus2.csv', delimiter=',', dtype=np.float64)
#stress2=np.loadtxt('OUT\yld_locus2.csv', delimiter=',', dtype=np.float64)
cmap_data = np.loadtxt('OUT\CMAP.csv', delimiter=',', dtype=np.float64)
debug = np.loadtxt('DEBUG\DEBUG.csv', delimiter=',', dtype=np.float64)

indx=0
#ang1=361
ang1=721
#ang2=181
ang2=361
cmap_output=np.zeros([ang1,ang2,3])
cmap_fail=np.zeros([ang1,ang2,3])
#########################################################################################
# PLOT PI-PLANE & WRITE OUTPUT DATA
# SET PLOTTING PARAMETERS #
print_par=1       # Print the result on csv file. (0: off / 1: on)
#r=np.max(stress[:,0])*1.5             # RADIUS OF AXIS
r=1.5
# OUTLINE OF PI-PLANE
# [1] PLOT YIELD_LOCUS
th=0*np.arccos(-1)/180
yl_axis=np.zeros([2,2])
x=r/np.sqrt(1+np.tan(th)**2)
yl_axis[0,0]=r/np.sqrt(1+np.tan(th)**2)
yl_axis[0,1]=np.tan(th)*yl_axis[0,0]
th=90*np.arccos(-1)/180
yl_axis[1,0]=r/np.sqrt(1+np.tan(th)**2)
yl_axis[1,1]=np.tan(th)*yl_axis[1,0]
r2=r

plt.figure(1, figsize=(8,8))
plt.grid(which='major', alpha=0.5)
plt.axis([-r2,r2,-r2,r2])
#plt.arrow(-yl_axis[0,0],yl_axis[0,1],2*yl_axis[0,0],2*yl_axis[0,1], width=0.006, color='black')
#plt.arrow(yl_axis[1,0],-yl_axis[1,1],2*yl_axis[1,0],2*yl_axis[1,1], width=0.006, color='black')
#plt.annotate("$\sigma_x$", xy=(yl_axis[0,0],yl_axis[0,1]), size=20, xytext=(yl_axis[0,0]-0.2,yl_axis[0,1]-0.15))
#plt.annotate("$\sigma_y$", xy=(yl_axis[1,0],yl_axis[1,1]), size=20, xytext=(yl_axis[1,0]-0.2,yl_axis[1,1]-0.15))
plt.axvline(x=0, color='k', linewidth=1.5)
plt.axhline(y=0, color='k', linewidth=1.5)
plt.tick_params(axis='both', direction='in', length=5, pad=6, labelsize=14)

# Original yield locus
plt.plot(stress[:,2], stress[:,3], 'k--')
# HAH20h
plt.plot(stress[:,0], stress[:,1], linewidth=2.5, c='black')
# Biaxial tension
len0=0
pv=2
px=(len0+pv)*0.5
py=-px+pv
#plt.plot([px,py], [py,px], linewidth=2.5)
# Plane strain tension
len0=0.7
pv=1.6
px=(len0+pv)*0.5
py=-px+pv
plt.plot([px,py], [py,px], linewidth=2.5)
# Uniaxial tension
len0=1
pv=1
px=(len0+pv)*0.5
py=-px+pv
plt.plot([px,py], [py,px], linewidth=2.5)
# Pure shear
len0=1.15
pv=0
px=(len0+pv)*0.5
py=-px+pv
plt.plot([px,py], [py,px], linewidth=2.5)
# Uniaxial compression
len0=1
pv=-1
px=(len0+pv)*0.5
py=-px+pv
plt.plot([px,py], [py,px], linewidth=2.5)
# Plane strain compression
len0=0.7
pv=-1.6
px=(len0+pv)*0.5
py=-px+pv
plt.plot([px,py], [py,px], linewidth=2.5)
# Balanced biaxial compression
pv=-2
px=(len0+pv)*0.5
py=-px+pv
#plt.plot([px,py], [py,px], linewidth=2.5)
plt.xlabel('$\sigma_x$/$\sigma_0$', size=20)
plt.ylabel('$\sigma_y$/$\sigma_0$', size=20)
#plt.legend()
plt.show()


