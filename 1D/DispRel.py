#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Illustrative justification of the use of the Bloch formalism in periodic media
It just plots the dispersion relation in 1D under different points of views.
"""
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Circle
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = 18
width = 2
plt.close("all")

#path_save="/home/philippe/Documents/ESA/Python_SC/1D/Pics/"

c = 343.2
pnts = 1500
k = np.linspace(-10,10,pnts)
k1 = np.linspace(-np.pi,np.pi,pnts)
k2m = np.linspace(-2*np.pi,-np.pi,pnts)
k2p = np.linspace(np.pi,2*np.pi,pnts)
k3m = np.linspace(-3*np.pi,-2*np.pi,pnts)
k3p = np.linspace(2*np.pi,3*np.pi,pnts)
k4m = np.linspace(-3*np.pi,-4*np.pi,pnts)
k4p = np.linspace(3*np.pi,4*np.pi,pnts)

omega = k*c
dispRel = np.abs(omega)
fig,ax = plt.subplots(1,figsize=(9, 4))

ax.plot(k,dispRel,"k",linewidth = width)
#ax.plot(k1,np.abs(k1*c),"r",linewidth = width)
#ax.plot(k2m,np.abs(k2m*c),"g",linewidth = width)
#ax.plot(k2p,np.abs(k2p*c),"g",linewidth = width)
#ax.plot(k3m,np.abs(k3m*c),"b",linewidth = width)
#ax.plot(k3p,np.abs(k3p*c),"b",linewidth = width)
##
#I = [-3,-2,-1,1,2,3]
#for i in I:
#    ax.plot(i*np.pi*np.ones(len(k)),omega,"k--",alpha=0.5)
ax.set_xlabel("$k$",fontsize=font)

ax.set_ylabel("$\omega$",fontsize=font,rotation=0)
ax.yaxis.set_label_coords(0.5, 1.02)
ax.grid(True,alpha=0.0)
    

ax.tick_params(labelsize=font)
ax.set_ylim([0,3000])
ax.set_xlim([-9,9])
    
ax.spines['left'].set_position('center')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
#plt.savefig(path_save+"disp1.pdf", bbox_inches='tight')

fig,ax = plt.subplots(1,figsize=(3, 4))
#
ax.plot(k1,np.abs(k1*c),"r",linewidth = width)
ax.plot(k1,-(np.abs(k1*c)-2*np.pi*c),"g",linewidth = width)
ax.plot(k1,np.abs(k1*c)+2*np.pi*c,"b",linewidth = width)
ax.plot(-np.pi*np.ones(len(k)),omega,"k--",alpha=0.5)
ax.plot(np.pi*np.ones(len(k)),omega,"k--",alpha=0.5)
#
ax.set_ylabel("$\omega$",fontsize=font-2,rotation=0)
ax.yaxis.set_label_coords(0.5, 1.02)
ax.set_xlabel("$k$",fontsize=font-2)
ax.grid(True,alpha=0.0)
#    
xT = [-np.pi,-np.pi/2,0,np.pi/2,np.pi]
labels = ['$-\pi$','$-\pi/2$', '0','$\pi/2$', '$\pi$']
plt.xticks(xT, labels,fontsize=font-2)
#
ax.tick_params(labelsize=font-2)
#
ax.set_ylim([0,3000])
ax.set_xlim([-np.pi,np.pi])

#
ax.spines['left'].set_position('center')
#
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
#
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.show()
#plt.savefig(path_save+"dispBloch.pdf", bbox_inches='tight')