#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute the plane wave expansion for a 2d sonic crystal. 
Code slightly modified from PhD thesis of Daniel Elford: Band gap formation in acoustically resonant phononic crystals.
Plot the band gap structure

The code only works for the circle case with a lattice constant equal to 1...
"""
import numpy as np
import scipy.special as sp
from numpy.linalg import inv
from scipy.linalg import eig
import matplotlib.pyplot as plt

# plotting parameters
plt.rc('text', usetex=True)
plt.rc('font', family='serif')  
font = 14
save_figs = 0


num_Eigenvals = 10 # number of eigenvalues to be computed
nmax = 7 # number of coefficients in the fourier expansion 

n1=np.arange(-nmax,nmax+1) 
n2=np.arange(-np.floor((2*nmax+1)**2/2),np.floor((2*nmax+1)**2/2)+1)

lat_a=1 # lattice constant
ro=0.4 # circle radius
l = 0.7 # side of the square for square inclusions
case = "circle" # type of inclusion
name = "circle"+str(int(100*ro))
# define the "IBZ" for PWE
K1=[np.array([i/lat_a,0.]) for i in np.arange(0.01,0.51,0.01)] 
K2=[np.array([0.5,i/lat_a]) for i in np.arange(0.01,0.51,0.01)]
K3=[np.array([i/lat_a,i/lat_a]) for i in np.arange(0.01,0.51,0.01)[::-1]]
K=K1+K2+K3    
         
# paramters of the host medium 
ch = 343.2 # speed of sound host medium
rhoh = 1.204  # density host medium
Bh=(ch**2)*rhoh
# paramters of the inclusion
rhoa=7800
ca=6100
Ba=(ca**2)*rhoa

#compute acoustic reflection coefficient
R=((ca*rhoa/(ch*rhoh)) - 1.0)/((ca*rhoa/(ch*rhoh)) + 1.0)
print("reflection coefficient={0}".format(R))

# generate grid
nx_arr=np.tile(n1,(1,2*nmax+1))
ny_arr=np.round(n2/(2*nmax+1))
[Nx,Nxp]=np.meshgrid(nx_arr, nx_arr)
[Ny,Nyp]=np.meshgrid(ny_arr, ny_arr)

fill = np.ones([1,np.size(nx_arr)])
Id_N = np.identity(np.size(nx_arr))

if case == "circle":
    f = np.pi*ro**2/lat_a**2
    #G_Gp = np.sqrt(4*np.pi*f*((Nx-Nxp)**2+(Ny-Nyp)**2)) + Id_N
    G_Gp=2*np.pi*ro/lat_a*np.sqrt((Nx-Nxp)**2+(Ny-Nyp)**2) + Id_N
    fG_Gp=2*f*sp.jvp(1,G_Gp,0)/G_Gp
elif case == "square": #not working...
    f = l**2 / lat_a**2
    Gx_p = 2*np.pi*l*(Nx-Nxp) + Id_N
    Gy_p = 2*np.pi*l*(Ny-Nyp) + Id_N #2*np.pi/lat_a
    fG_Gp = f*np.sinc(Gx_p*l/2)*np.sinc(Gy_p*l/2) #(Nx-Nxp), Nx
else:
    raise ValueError('not implemented yet !')
 
delp=(rhoh/rhoa-1)/(f*rhoh/rhoa+1-f) 
delt=(Bh/Ba-1)/(f*Bh/Ba+1-f)
    

BG = np.zeros([num_Eigenvals,len(K)],dtype=complex)
print("Main loop begin")               
for i in range(len(K)):
    kx,ky = K[i]

    mag_k_p_g=((kx+Nx)**2+(ky+Ny)**2) #Gamma_1 + Gamma_2
    k_p_g_dot_kp_p_g= (kx+Nx)*(kx+Nxp)  +  (ky+Ny)*(ky+Nyp)
                  
    kro_del_ggp=(Nxp==Nx)&(Nyp==Ny)
            
    M = mag_k_p_g * kro_del_ggp + delp * fG_Gp * k_p_g_dot_kp_p_g * (1-kro_del_ggp)
    N =  kro_del_ggp + fG_Gp * delt * (1-kro_del_ggp)
            
    A=np.dot(inv(N),M)
    b,a=eig(A)
    g=b[b>0]
    eigs=np.sort(g)
            
    BG[:,i] = eigs[0:num_Eigenvals]
    k=np.tile(kx,(1,num_Eigenvals)) #kxt
    min(b)
    if i%5 == 0:
        print("Progress ",i,"/",len(K))

################## plot results ###################################            
plt.figure(num=5,figsize=(4,6))
for i in range(0,10):
    ydata = ch*np.real(np.sqrt(BG[i,:]))
    #plt.scatter(range(ydata),ydata,marker='+',linewidths =0.5,color="k",alpha=1)
    plt.plot(ydata,linewidth=2,color="k",alpha=0.25)
    
plt.grid(True)
plt.xticks([0,len(K1) ,len(K1)+len(K2) ,len(K1)+len(K2)+len(K3)],["$\Gamma$","$X$","$M$","$\Gamma$"],fontsize=font)
plt.ylabel("Frequency [Hz]",fontsize=font+2)
plt.xlabel("Bloch Wavevector",fontsize=font+2)
plt.yticks(fontsize=font+2)
plt.ylim(0,500)
plt.xlim(0,len(K1)+len(K2)+len(K3))  
plt.show()

if save_figs ==1:
    plt.savefig('/home/philippe/Documents/ESA/Report/fig/2D/BG/mesh/2D_PWEBG_'+ name +'bis.pdf', bbox_inches='tight')