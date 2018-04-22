#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 11:56:42 2018

@author: philippe
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from scipy.integrate import simps
import matplotlib.colors as colors
#from matplotlib2tikz import save as tikz_save


def TMM(rho,c,h,k,theta): # construct one layer of the unit cell
    T11 = np.cos(k*h*np.cos(theta))
    T12 = 1j*rho*c*np.sin(k*h*np.cos(theta)) / np.cos(theta)
    T21 = 1j*np.cos(theta)*np.sin(k*h*np.cos(theta)) / (rho*c)
    T22 = np.cos(k*h*np.cos(theta))
    T = np.array([[T11, T12],[T21, T22]])
    return T

def MatUnitCell3(Ta,Tb,Tc,freq): # 3 layers unit cell
    Tu = np.zeros([2,2,len(freq)], dtype=np.complex)
    for i in range(0,len(freq)-1):
        Tu[:,:,i] = np.dot(np.dot(Ta[:,:,i],Tb[:,:,i]),Tc[:,:,i])   
        
    Zu = Tu[0,0,:]/Tu[1,0,:]
    keff = (1/hu)*np.arccos(0.5*(Tu[0,0,:]+Tu[1,1,:]))
    return Tu,Zu,keff

def MatUnitCell2(Ta,Tb,freq): # 2 layers unit cell
    Tu = np.zeros([2,2,len(freq)], dtype=np.complex)
    for i in range(0,len(freq)-1):
        Tu[:,:,i] = np.dot(Ta[:,:,i],Tb[:,:,i])   
        
    Zu = Tu[0,0,:]/Tu[1,0,:]
    keff = (1/hu)*np.arccos(0.5*(Tu[0,0,:]+Tu[1,1,:]))
    return Tu,Zu,keff

def MultiLayerTMM(Tu,Nlayers):
    k = 2
    Tg = np.zeros([2,2,len(freq)], dtype=np.complex)
    
    if Nlayers == 1:
        Tg = Tu
        return Tg
    
    for i in range(0,len(freq)-1):
        Tg[:,:,i] = np.dot(Tu[:,:,i],Tu[:,:,i])   
       
    TgMult = Tg
    while k < Nlayers:
        for i in range(0,len(freq)-1):
            TgMult[:,:,i] = np.dot(Tg[:,:,i],Tu[:,:,i])   
        k = k + 1
    Tg = TgMult     
    return Tg

def ReflectionCoeff(theta,Zhost,TG):
    R = (TG[0,0,:] + TG[0,1,:]/(np.cos(theta)*Zhost) - Zhost*TG[1,0,:]*np.cos(theta) - TG[1,1,:]) \
     / ( - TG[0,0,:] + TG[0,1,:]/(np.cos(theta)*Zhost) + Zhost*TG[1,0,:]*np.cos(theta) - TG[1,1,:])     
    R = np.abs(R) 
    return R

def TransmissionLoss(theta,Zhost,TG): # sur ou fois cos theta ?
    
    T = np.abs(2/( -TG[0,0,:] - TG[1,1,:] + Zhost*TG[1,0,:]*np.cos(theta) + TG[0,1,:]/(np.cos(theta)*Zhost) )    )
    TL = np.abs(TG[0,0,:] + TG[0,1,:]*np.cos(theta)/Zhost + Zhost*TG[1,0,:]/np.cos(theta) + TG[1,1,:])**2/4
    TLdB = 10*np.log10(TL)
    return TLdB, TL, T

def DiffuseTL(Nlayers):
    theta = np.linspace(0.01,np.pi/2,len(freq))
    TLDif = np.zeros([len(freq),len(freq)])
    TLDif_red = np.zeros([len(freq)])
    Tu_tmp = np.zeros([2,2,len(freq)])
    Ttheta = np.zeros([2,2,len(freq)])
    
    for i in range(0,len(theta)-1):
        Tu_tmp = MatUnitCell3(TMM(rhoa,ca,ha,ka,theta[i]),TMM(rhow,cw,hw,kw,theta[i]),TMM(rhoa,ca,ha,ka,theta[i]),freq)[0]
        Ttheta = MultiLayerTMM(Tu_tmp,Nlayers)
        TLDif[:,i] = TransmissionLoss(theta[i],Zf,Ttheta)[0]*np.cos(theta[i])*np.sin(theta[i])
    
    for j in range(0,len(theta)-1):
        TLDif_red[j] = 2*sc.integrate.simps(TLDif[j,:],theta)
        
    return TLDif_red
################## Definition of the frequency domain ######################
N_samp = 400 # number of sampling points for frequency and incidence angle
freq_ini = 1
freq_end = 8e3
freq = np.linspace(freq_ini,freq_end,N_samp)
omega = 2*np.pi*freq
############################################################

########### angle of incidence of the plane wave ###########
theta = 0
############################################################

########## Parameters of the host medium ##################
cf = 343.2 #1500. # speed of sound of the host medium
rhof = 1.204 #1000.  # density of the host medium
Zf = cf*rhof
kf = 2*np.pi*freq/cf
############################################################

############ definition of the unit cell ####################
# layer 1 , acier
Ea = 200E9 #69
rhoa = 7800 #2700
ca = np.sqrt(Ea/rhoa)*(1+1j*0.)#6200
ha = 1e-2 #5e-3
ka = 2*np.pi*freq/ca
###################################
# layer 2 
cb = 1000*(1+1j*0.03)
rhob = 1250
hb = 30e-2
kb = 2*np.pi*freq/cb
###################################
# helmholtz resonator layer - Porous model

#c_res = 
#rho_res =
#h_res = 10e-2
###################################
#water layer
cw = 1481*(1+1j*0.03)
cw0 = 1481*(1+1j*0.)
rhow = 1000
hw = 25e-2
Zw = cw*rhow
kw = 2*np.pi*freq/cw
kw0 = 2*np.pi*freq/cw0
###################################
# porous layer - Darcy Law
Keq = 1.4*101325
sigma = 100
rhoeq = rhof*(1+sigma/(1j*omega*rhof))
ceq = np.sqrt(Keq/rhoeq)
Zeq = rhoeq*ceq
keq = 2*np.pi*freq/ceq
heq = 0.005e-3

#plt.figure(1)
#plt.subplot(211)
#plt.plot(freq, np.abs(rhoeq), 'b--')
#plt.grid(True)
#plt.subplot(212)
#plt.plot(freq, np.abs(ceq), 'r--')
#plt.grid(True)
#plt.show()


############################################################
hu = (ha + hw) # global thickness of the multi-layer

Ta = TMM(rhoa,ca,ha,ka,theta)
Tb = TMM(rhob,cb,hb,kb,theta)
Tw = TMM(rhow,cw,hw,kw,theta)
Tw0 = TMM(rhow,cw0,hw,kw0,theta)

Tu,Zu,keff = MatUnitCell3(Ta,Tw,Ta,freq)
Tu0,Zu0,keff0 = MatUnitCell3(Ta,Tw0,Ta,freq)
#Tu,Zu,keff = MatUnitCell3(Ta,Teq,Ta,freq)


ceff = 2*np.pi*freq / keff
rhoeff = 1j*np.sin(keff*hu) / (ceff*Tu[1,0])
#rhoeff = (ha*rhoa + hb*rhob)/(hu)
Zeff = rhoeff*ceff
#Teff = TMM(rhoeff,ceff,ha+hb,keff,theta)

TG = MultiLayerTMM(Tu,1)
TG3 = MultiLayerTMM(Tu,3)
TG10 = MultiLayerTMM(Tu,10)

Zeff1 = TG[0,0,:]/TG[1,0,:]
Zeff3 = TG3[0,0,:]/TG3[1,0,:]
Zeff10 = TG10[0,0,:]/TG10[1,0,:]

TL = TransmissionLoss(theta,Zf,TG)[0]
T = TransmissionLoss(theta,Zf,TG)[2]
R = ReflectionCoeff(theta,Zf,TG)
A = 1 - T**2 - R**2

TL3 = TransmissionLoss(theta,Zf,TG3)[0]
TL10 = TransmissionLoss(theta,Zf,TG10)[0]
#TLeff = TransmissionLoss(theta,Zf,Teff)

TLDif = DiffuseTL(1)
TL3Dif = DiffuseTL(3)
TL10Dif = DiffuseTL(10)

######################### plot figures ###############################
a = 1.5
fontsize = 12
plt.figure(num=2, figsize=(6,8))

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#plt.rc('font', **{'family': 'serif', 'serif': ['impact']}) #Computer Modern

plt.subplot(221)
plt.plot(np.real(keff*hu),freq,color='red',label='3 \% damping - water',linewidth=1.2*a)
plt.plot(np.real(keff0*hu),freq,color='darkred',linestyle='--', label='no damping',linewidth=0.8*a)
plt.ylabel('Frequency (Hz)',fontsize=fontsize)
plt.xlabel('$ ka $',fontsize=fontsize)
plt.axhspan(1905.0,2950.0,alpha=0.2, color = 'grey')
plt.axhspan(4153.0,5911.0,alpha=0.2, color = 'grey')
plt.axhspan(6736.0,7900.0,alpha=0.2, color = 'grey')
plt.xlim(0,np.pi*0.99)
plt.ylim(0,7900)
plt.xticks(fontsize=fontsize-2)
plt.yticks(fontsize=fontsize-2)
plt.grid(True)
plt.legend(loc='best',fontsize=fontsize-4)

plt.subplot(222)
plt.plot(TL,freq,color='blue',label='L=1',linewidth=0.9*a, linestyle='--')
plt.plot(TL3,freq,color='steelblue',label='L=3',linewidth=0.8*a, linestyle='-.')
plt.plot(TL10,freq,color='navy',label='L=10',linewidth=0.7*a, linestyle=':')
#plt.plot(TL15,freq,color='black',label='L=15',linewidth=a)
#plt.ylabel('Frequency (Hz)',fontsize=font)
plt.xlabel('TL (dB)',fontsize=fontsize)
plt.xlim(45,130)
plt.ylim(0,7900)
plt.axhspan(1905.0,2950.0,alpha=0.2, color = 'grey')
plt.axhspan(4153.0,5911.0,alpha=0.2, color = 'grey')
plt.axhspan(6736.0,7900.0,alpha=0.2, color = 'grey')
plt.xticks(fontsize=fontsize-2)
plt.yticks(fontsize=fontsize-2)
plt.grid(True)
#plt.legend(loc='best')

plt.subplot(224)
plt.plot(TLDif,freq,color='blue',label='L=1',linewidth=0.9*a, linestyle='--')
plt.plot(TL3Dif,freq,color='steelblue',label='L=3',linewidth=0.8*a,linestyle='-.')
plt.plot(TL10Dif,freq,color='navy',label='L=10',linewidth=0.7*a, linestyle=':')
#plt.plot(TL15Dif,freq*1e-3,color='black',label='L=15',linewidth=a)
#plt.ylabel('Frequency (Hz)',fontsize=font)
plt.xlabel('TL diffuse field (dB)',fontsize=fontsize)
plt.xlim(45,130)
plt.ylim(0,7900)
plt.axhspan(1905.0,2950.0,alpha=0.2, color = 'grey')
plt.axhspan(4153.0,5911.0,alpha=0.2, color = 'grey')
plt.axhspan(6736.0,7900.0,alpha=0.2, color = 'grey')
plt.xticks(fontsize=fontsize-2)
plt.yticks(fontsize=fontsize-2)
plt.grid(True)
#plt.yticks([])
#plt.legend(loc='best')
#plt.savefig('TLdiffuse.png', bbox_inches='tight')

plt.subplot(223)

plt.plot(np.abs(Zeff/Zf),freq,color='red',label='$L=\infty$',linewidth=1.2*a)
plt.plot(np.abs(Zeff1/Zf),freq,color='blue',label='$L=1$',linewidth=0.9*a, linestyle='--')
plt.plot(np.abs(Zeff3/Zf),freq,color='steelblue',label='$L=3$',linewidth=0.8*a, linestyle='-.')
plt.plot(np.abs(Zeff10/Zf),freq,color='navy',label='$L=10$',linewidth=0.7*a, linestyle=':')
plt.ylabel('Frequency (Hz)',fontsize=fontsize)
plt.xlabel('$Z_{eff}/Z_{air}$',fontsize=fontsize)
plt.ylim(0,7900)
plt.xlim(0,30000)
plt.axhspan(1905.0,2950.0,alpha=0.2, color = 'grey')
plt.axhspan(4153.0,5911.0,alpha=0.2, color = 'grey')
plt.axhspan(6736.0,7900.0,alpha=0.2, color = 'grey')
plt.xticks(fontsize=fontsize-2)
plt.yticks(fontsize=fontsize-2)
plt.grid(True)
plt.legend(loc=1,fontsize=fontsize-4)

plt.savefig('/home/philippe/Documents/ESA/Python_SC/1D/Pics/1D_media3.pdf', bbox_inches='tight')
#plt.savefig('/home/philippe/Documents/ESA/Report/fig/1D/1D_media3.pdf', bbox_inches='tight')
#plt.savefig('/home/philippe/Documents/ESA/Report/fig/1D/1D_media.eps', bbox_inches='tight')
#plt.savefig('/home/philippe/Documents/ESA/Report/fig/1D/1D_media.png', bbox_inches='tight')

#tikz_save()
#tikz_save('mytikz.tex',figureheight='\\figureheight',figurewidth='\\figurewidth')