# -*- coding: utf-8 -*-
"""
Sound propagation in periodic media - 1D analytic cases

2 small examples:
    mass_spring: two masses of weight m1 and m2 connected by a spring. This pattern is repeted to created an infinite 1D periodic medium.
    The plot shows the 2 eigenvalues of the problem - as a function of the wavenumber k. 
    The slope of the dispersion relation represents the speed of the wave propagation.
    
    Helmholtz: acoustic absorption coefficient of a Helmholtz resonator as a fonction of the air temperature. The model is taken from this paper:
        https://doi.org/10.1016/j.jsv.2007.01.012
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
plt.close('all')
font = 11
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


# set the below variables to True to see the related plots
mass_spring = True
Helmholtz = True
######################

if mass_spring == True:
    a = 3 # lattice parameter
    beta = 2 # spring constant
    m1 = 5 # 1st mass
    m2 = 3 #2nd mass
    k = np.linspace(0,np.pi/(2*a),500)
    
    omega1 = np.sqrt(   beta*(1/m1 + 1/m2) + np.sqrt( (beta*(1/m1 + 1/m2))**2  - (2*beta*np.sin(k*a))**2 / (m1*m2) ) )
    omega2 = np.sqrt(   beta*(1/m1 + 1/m2) - np.sqrt( (beta*(1/m1 + 1/m2))**2  - (2*beta*np.sin(k*a))**2 / (m1*m2) ) )
    

    plt.figure(1)
    plt.plot(k, omega1, linestyle='--', color='r', label='$\omega_1$')
    plt.plot(k, omega2, linestyle='--', color='b', label='$\omega_2$')
    plt.xlabel('$k$',fontsize = font)
    plt.ylabel('$\omega(k)$',fontsize = font)
    plt.legend(fontsize = font)
    plt.xticks(fontsize = font)
    plt.yticks(fontsize = font)
    plt.grid(True)
    plt.xlim(0,max(k))
    plt.ylim(0,1.1*max(max(omega1),max(omega2)))
    plt.show()


if Helmholtz == True:
    
    # set the frequency range
    min_freq = 1
    max_freq = 15000
    N_freq = max_freq
    freq = np.linspace(min_freq,max_freq,N_freq)
    omega = 2*np.pi*freq
    
    # define parameters
    kappa = 24.36e-3 # thermal conductivity 0 degree
    Cp = 1.006 # isobar heat capacity 
    gamma = 1.4 # heat capacity ratio
    Patm = 101325 # atmospheric pressure
    T0 = 273.15
    c0 = 331.4 # sound speed air 0 degree
    eta0 = 1.71e-5 # air dynamic viscosity 0 degree
    C = 120 # Sutherland's constant
    rho0 = 1.292 # air density 0 degree
    Z0 = rho0*c0

    T = T0 + np.arange(0,600,100)
    N_T = len(T)
    c = 20.05*np.sqrt(T)
    eta = eta0*((T0+C)/(T+C))*(T/T0)**(1.5)
    rho = rho0*T0/T
    Z = rho*c

    # resonator parameters, in meters
    d = 0.8e-3# length of the neck
    phi = 0.05# porosity or open area ratio
    r = 0.15e-3# radius of the neck
    eps_0 = (8*r)/(3*np.pi)
    eps_e = eps_0*(1 - 1.14*np.sqrt(phi))
    L = 0.02#size of the cavity

    # equivalent porous parameters for T[0]
    sigma = 8*eta/(phi*r**2)
    tortuosity = 1 + 2*eps_e/d
    
    GJ = np.sqrt(1+1j*(4*omega*rho[0]*tortuosity**2*eta[0])/((sigma[0]*phi*r)**2))
    rho_eq = (1/phi)*rho[0]*tortuosity*(1+sigma[0]*phi*GJ[0]/(1j*omega*rho[0]*tortuosity)) # equivalent density
    
    GJ2 = np.sqrt(1+1j*omega*(r**2*Cp*rho[0]/(16*kappa)))
    beta = gamma-(gamma-1)*(1+8*kappa*GJ2[0]/(1j*omega*r**2*Cp*rho[0])) # equivalent compressibility
    
    K_eq = gamma*Patm/(phi*beta)
    c_eq = np.sqrt(K_eq/rho_eq) # effective speed of sound at T[0]
    
    # intialize variables
    k = np.zeros((N_T,N_freq))
    Rs = np.zeros((N_T,N_freq))
    ZB = np.zeros((N_T,N_freq),dtype=np.complex)
    R = np.zeros((N_T,N_freq),dtype=np.complex)
    alpha = np.zeros((N_T,N_freq))
    
    # loop over the different temperatures
    for t in range(N_T):
        k[t,:] = omega/c[t]
        Rs[t,:] = 0.5*np.sqrt(2*eta[t]*omega*rho[t]) #surface resistance
        ZB[t,:] = (2*d + 4*eps_e)*(Rs[t,:]/(phi*r)) + 1j*(  (1/phi)*omega*rho[t]*(2*eps_e + d) 
            - rho[t]*c[t]*(1/np.tan(k[t,:]*L)) ) # impedance at the plate surface
        
        R[t,:] = (ZB[t,:] - Z[t]) / (ZB[t,:] + Z[t]) # reflection coefficient
        alpha[t,:] = 1 - abs(R[t,:])**2 # absorption coefficient
        
        
    # estimate the resonance frequencies
    freq_res = 1/(2*np.pi) * np.sqrt(  (c**2*phi)/  ((2*eps_e+d)*L) ) 
    print("Air temperature:", T[0], "K, approximate resonance frequency = %.1f" %freq_res[0], "Hz")

    # plot the absorption coefficient for different temperatures
    NUM_COLORS = N_T
    cm = plt.get_cmap('jet')
    fig = plt.figure(figsize=(12,4))
    ax = fig.add_subplot(111)
    ax.set_prop_cycle(color=[cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
    for i in range(NUM_COLORS):
        ax.plot(freq,alpha[i,:],linestyle='-',linewidth=2.0)
        plt.ylabel(r"$\alpha$",fontsize = font)
        plt.xlabel('Frequency (Hz)',fontsize = font)
        plt.xlim(0,max_freq)
        plt.ylim(0,1)    
        #plt.xticks([0.0,np.pi/2,np.pi],["$0$","$\pi/2$","$\pi$"],fontsize = font)
        plt.xticks(fontsize = font)    
        plt.yticks(fontsize = font)
        plt.grid(True)
        
        
    norm = matplotlib.colors.Normalize(vmin=np.min(T),vmax=np.max(T))
    s_m = matplotlib.cm.ScalarMappable(cmap=cm,norm=norm)
    s_m.set_array([])
    p = fig.colorbar(s_m)
    p.set_label('Temperature (K)', rotation=90, fontsize = font)
    plt.show()