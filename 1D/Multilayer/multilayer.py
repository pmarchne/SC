# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 16:35:48 2018

@author: Philippe Marchner

Band Gap computation of a 1D periodic media with losses
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib
import scipy as sc
from scipy.integrate import simps
plt.close("all")
################## Frequency domain ######################
N_samp = 500 # number of sampling points
freq_ini = 1
freq_end = 8e3
freq = np.linspace(freq_ini,freq_end,N_samp)
omega=2*np.pi*freq


########## Parameters of the host medium, air ##################
cf = 343.2 # speed of sound of the host medium
rhof = 1.204  # density of the host medium
Zf = cf*rhof # characteristic impedance
kf = 2*np.pi*freq/cf # wavenumber
theta = 0 # incidence angle of the plane wave
############################################################


############ definition of the unit cell, composed of two sublayers ####################
# layer 1, steel
Ea = 200E9
rhoa = 7800
ca = np.sqrt(Ea/rhoa) 

Za = ca*rhoa
ha = 1e-2#15e-3 # width of the first layer
ka = 2*np.pi*freq/ca
# define the tranfer matrix
T11a = np.cos(ka*ha*np.cos(theta))
T12a = 1j*rhoa*ca*np.sin(ka*ha*np.cos(theta)) / np.cos(theta)
T21a = 1j*np.cos(theta)*np.sin(ka*ha*np.cos(theta)) / (rhoa*ca)
T22a = np.cos(ka*ha*np.cos(theta))
Ta = np.array([[T11a, T12a],[T21a, T22a]])

# layer 2, water
cb = 1481#1000*(1+1j*0.00) #0.02
rhob = 1000#1250
Zb = cb*rhob
hb = 25e-2#15e-3 #30e-3
kb = 2*np.pi*freq/cb
# define the tranfer matrix
T11b = np.cos(kb*hb*np.cos(theta))
T12b = 1j*rhob*cb*np.sin(kb*hb*np.cos(theta)) / np.cos(theta)
T21b = 1j*np.cos(theta)*np.sin(kb*hb*np.cos(theta)) / (rhob*cb)
T22b = np.cos(kb*hb*np.cos(theta))
Tb = np.array([[T11b, T12b],[T21b, T22b]])
############################################################
hu2 = ha + hb # global thickness of the unit cell
hu3 = 2*ha + hb

# tranfer matrix of the unit cell, 3 layers
Tu3 = np.zeros([2,2,len(freq)], dtype=np.complex)
for i in range(0,len(freq)-1):
    Tu3[:,:,i] = np.dot(np.dot(Ta[:,:,i],Tb[:,:,i]),Ta[:,:,i])
    
# tranfer matrix of the unit cell, 3 layers
Tu2 = np.zeros([2,2,len(freq)], dtype=np.complex)
for i in range(0,len(freq)-1):
    Tu2[:,:,i] = np.dot(Ta[:,:,i],Tb[:,:,i])


# compute effective wavenumber 
keff2 = (1/hu2)*np.arccos(0.5*(Tu2[0,0,:]+Tu2[1,1,:]))
keff3 = (1/hu3)*np.arccos(0.5*(Tu3[0,0,:]+Tu3[1,1,:]))

am = np.array([0.0,0.02,0.04,0.06,0.08,0.1])
cb = 1481*(1+1j*am)
Nrho = len(cb)
gamma = np.zeros(Nrho,dtype=np.complex)
keff_am = np.zeros([Nrho,len(freq)],dtype=np.complex)
for i in range(Nrho):
    gamma[i] = (rhoa*ca)/(rhob*cb[i])
    keff_am[i,:] = (1/hu2)*np.arccos(np.cos(omega*ha/ca)*np.cos(omega*hb/cb[i]) - 0.5*(1/gamma[i] + gamma[i])*np.sin(omega*ha/ca)*np.sin(omega*hb/cb[i]))


ceff2 = 2*np.pi*freq / keff2
rhoeff2 = 1j*np.sin(keff2*hu2) / (ceff2*Tu2[1,0])
Zeff2 = rhoeff2*ceff2

ceff3 = 2*np.pi*freq / keff3
rhoeff3 = 1j*np.sin(keff3*hu3) / (ceff3*Tu3[1,0])
Zeff3 = rhoeff3*ceff3

######################### plot figures ###############################
font = 11
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.figure(2, figsize=(5,5))
plt.subplot(121)
plt.plot(np.real(keff2*hu2),freq*1e-3,color="red",label="$h_1$ = 1 cm")
plt.plot(np.real(keff3*hu3),freq*1e-3,color="red",linestyle='--',label="$h_1$ = 2 cm")
plt.ylabel('Frequency [kHz]',fontsize = font)
plt.xlabel('Re($ka$)',fontsize = font)
plt.xlim(0,np.pi*0.99)
plt.xticks([0.0,np.pi/2,np.pi],["$0$","$\pi/2$","$\pi$"],fontsize = font)
plt.ylim(0,8)
plt.xticks(fontsize = font)
plt.yticks(fontsize = font)
plt.grid(True)

plt.subplot(122)
plt.plot(np.abs(np.imag(keff2*hu2)),freq*1e-3,color="red",label="$h_1$ = 1 cm")
plt.plot(np.abs(np.imag(keff3*hu3)),freq*1e-3,color="red",linestyle='--',label="$h_1$ = 2 cm")
plt.xticks(fontsize = font)
#plt.xticks([0.0,np.pi/2,np.pi],["$0$","$\pi/2$","$\pi$"],fontsize = font)
plt.xticks([0.0,np.pi/4,np.pi/2],["$0$","$\pi/4$","$\pi/2$"],fontsize = font)
plt.yticks(fontsize = font)
plt.xlabel('Im($ka$)',fontsize = font)
plt.xlim(-0.001,np.pi/2)
plt.ylim(0,8)
plt.grid(True)
plt.legend()

#plt.savefig('keff_red_slides.pdf', bbox_inches='tight')

plt.figure(3,figsize=(6,5))
plt.plot(freq*1e-3,np.abs(Zeff2/Zf),color='red',label="$h_1$ = 1 cm")
plt.plot(freq*1e-3,np.abs(Zeff3/Zf),color='red',linestyle='--',label="$h_1$ = 2 cm")
#plt.plot(freq,np.abs(Zeff3/Zf),color='red',linestyle='--')
plt.yscale("log")
plt.xticks(fontsize = font)
plt.yticks(fontsize = font)
plt.xlabel('Frequency [kHz]',fontsize = font)
plt.ylabel('$Z_{eff}/Z_{air}$',fontsize = font)
plt.ylim(100,100000)
plt.xlim(0,8)
plt.grid(True)
plt.legend()
#plt.legend()
#plt.savefig('Zeff_slides.pdf', bbox_inches='tight')


NUM_COLORS = Nrho

cm = plt.get_cmap('autumn')
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(121)
#ax.set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
ax.set_prop_cycle(color=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
for i in range(NUM_COLORS):
    ax.plot(np.real(keff_am[i]*hu2),freq*1e-3)
plt.ylabel('Frequency (kHz)',fontsize = font)
plt.xlabel('Re($ka$)',fontsize = font)
plt.xlim(0,np.pi)
plt.ylim(0,8)    
plt.xticks([0.0,np.pi/2,np.pi],["$0$","$\pi/2$","$\pi$"],fontsize = font)
plt.yticks(fontsize = font)
plt.grid(True)
  
ax = fig.add_subplot(122)
#ax.set_prop_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
ax.set_prop_cycle(color=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
for i in range(NUM_COLORS):
    ax.plot(np.abs(np.imag(keff_am[i]*hu2)),freq*1e-3)
 
plt.xlabel('Im($ka$)',fontsize = font)
plt.xlim(0,np.pi/2)
plt.ylim(0,8)    
plt.xticks([0.0,np.pi/4,np.pi/2],["$0$","$\pi/4$","$\pi/2$"],fontsize = font)
plt.yticks(fontsize = font)
plt.grid(True)


norm = matplotlib.colors.Normalize(vmin=np.min(am)*100,vmax=np.max(am)*100)
s_m = matplotlib.cm.ScalarMappable(cmap=cm,norm=norm)
s_m.set_array([])
p = fig.colorbar(s_m)
p.set_label('Damping (\%)', rotation=90, fontsize = font)
#plt.savefig('damping_slides.pdf', bbox_inches='tight')
#print(Za/Zb)



