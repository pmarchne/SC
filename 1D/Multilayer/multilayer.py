# -*- coding: utf-8 -*-
"""
Band Gap computation of a 1D periodic media with losses

parameters:
    - angle of the incidence wave
    - layer properties
    - unit cell configuration: 2-layers steel-water and steel-water-steel sandwich
"""

import numpy as np
np.seterr(divide='ignore', invalid='ignore')
import matplotlib.pyplot as plt
import matplotlib
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


############ definition of the unit cell ####################
# layer 1, steel
Ea = 200E9
rhoa = 7800
ca = np.sqrt(Ea/rhoa) 

Za = ca*rhoa
ha = 1e-2# width  of the first layer
ka = 2*np.pi*freq/ca
# layer 2, water 
cb = 1481
rhob = 1000
Zb = cb*rhob
hb = 25e-2
kb = 2*np.pi*freq/cb

# define the tranfer matrix
def TMM(rho,c,h,k,theta):
    T11 = np.cos(k*h*np.cos(theta))
    T12 = 1j*rho*c*np.sin(k*h*np.cos(theta)) / np.cos(theta)
    T21 = 1j*np.cos(theta)*np.sin(k*h*np.cos(theta)) / (rho*c)
    T22 = np.cos(k*h*np.cos(theta))
    return np.array([[T11, T12],[T21, T22]])
    
Ta = TMM(rhoa,ca,ha,ka,theta)
Tb = TMM(rhob,cb,hb,kb,theta)

# global thickness of the unit cell
hu2 = ha + hb # 2-layers, steel-water
hu3 = 2*ha + hb # sandwich layer, steel-water-steel
# because the situation is symmetric, the sandwich layer is equivalent to the 2-layers with a doubled steel thickness

# tranfer matrix of the unit cell, sandwich layer
Tu3 = np.zeros([2,2,len(freq)], dtype=np.complex)
for i in range(0,len(freq)-1):
    Tu3[:,:,i] = np.dot(np.dot(Ta[:,:,i],Tb[:,:,i]),Ta[:,:,i])
    
# tranfer matrix of the unit cell, 2 layers
Tu2 = np.zeros([2,2,len(freq)], dtype=np.complex)
for i in range(0,len(freq)-1):
    Tu2[:,:,i] = np.dot(Ta[:,:,i],Tb[:,:,i])


# compute effective parameters for the 2 configurations
keff2 = (1/hu2)*np.arccos(0.5*(Tu2[0,0,:]+Tu2[1,1,:]))
ceff2 = 2*np.pi*freq / keff2
rhoeff2 = 1j*np.sin(keff2*hu2) / (ceff2*Tu2[1,0])
Zeff2 = rhoeff2*ceff2

keff3 = (1/hu3)*np.arccos(0.5*(Tu3[0,0,:]+Tu3[1,1,:]))
ceff3 = 2*np.pi*freq / keff3
rhoeff3 = 1j*np.sin(keff3*hu3) / (ceff3*Tu3[1,0])
Zeff3 = rhoeff3*ceff3

# influence of damping in water on the band gap structure
am = np.array([0.0,0.02,0.04,0.06,0.08,0.1]) # damping in % 
cb = 1481*(1+1j*am)
Nrho = len(cb)
gamma = np.zeros(Nrho,dtype=np.complex)
keff_am = np.zeros([Nrho,len(freq)],dtype=np.complex)
for i in range(Nrho):
    # use the analytic dispersion relation
    gamma[i] = (rhoa*ca)/(rhob*cb[i])
    keff_am[i,:] = (1/hu2)*np.arccos(np.cos(omega*ha/ca)*np.cos(omega*hb/cb[i]) - 0.5*(1/gamma[i] + gamma[i])*np.sin(omega*ha/ca)*np.sin(omega*hb/cb[i]))


######################### plot figures ###############################
font = 16
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

Color = "red"
plt.figure(1, figsize=(5,5))
plt.subplot(121)
plt.plot(np.real(keff2*hu2),freq*1e-3,color=Color,label="$h_1$ = 1 cm")
plt.plot(np.real(keff3*hu3),freq*1e-3,color=Color,linestyle='--',label="$h_1$ = 2 cm")
plt.ylabel('Frequency [kHz]',fontsize = font)
plt.xlabel('Re($ka$)',fontsize = font)
plt.xlim(0,np.pi*0.99)
plt.xticks([0.0,np.pi/2,np.pi*0.99],["$0$","$\pi/2$","$\pi$"],fontsize = font)
plt.ylim(0,8*0.998)
plt.xticks(fontsize = font)
plt.yticks(fontsize = font)
plt.grid(True)

plt.subplot(122)
plt.plot(np.abs(np.imag(keff2*hu2)),freq*1e-3,color=Color,label="$h_1$ = 1 cm")
plt.plot(np.abs(np.imag(keff3*hu3)),freq*1e-3,color=Color,linestyle='--',label="$h_1$ = 2 cm")
plt.xticks(fontsize = font)
plt.xticks([0.0,np.pi/4,np.pi/2],["$0$","$\pi/4$","$\pi/2$"],fontsize = font)
plt.yticks(fontsize = font)
plt.xlabel('Im($ka$)',fontsize = font)
plt.xlim(-0.001,np.pi/2)
plt.ylim(0,8*0.998)
plt.grid(True)
plt.legend(fontsize = font-2)

#plt.savefig('keff_red_slides.pdf', bbox_inches='tight')

plt.figure(2,figsize=(6,5))
plt.plot(freq*1e-3,np.abs(Zeff2/Zf),color=Color,label="$h_s$ = 1 cm")
plt.plot(freq*1e-3,np.abs(Zeff3/Zf),color=Color,linestyle='--',label="$h_s$ = 2 cm")
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

cm = plt.get_cmap('magma')
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(121)
ax.set_prop_cycle(color=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
for i in range(NUM_COLORS):
    ax.plot(np.real(keff_am[i]*hu2),freq*1e-3)
plt.ylabel('Frequency [kHz]',fontsize = font)
plt.xlabel('Re($ka$)',fontsize = font)
plt.xlim(0,np.pi)
plt.ylim(0,8)    
plt.xticks([0.0,np.pi/2,np.pi],["$0$","$\pi/2$","$\pi$"],fontsize = font)
plt.yticks(fontsize = font)
plt.grid(True)
  
ax = fig.add_subplot(122)
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

plt.show()


