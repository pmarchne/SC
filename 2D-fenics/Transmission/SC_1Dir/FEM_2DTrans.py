#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solves the Helmholtz equation for a transmission problem with periodic boundary conditions.
It simulates the propagation of acoustic waves in a 1-direction finite sonic crystal (here the Gamma-X direction)

The code is only valid for low frequencies.
It computes the pressure maps for different frequencies and the transmission/reflection coefficients over a user-defined frequency range.

you can change the geometry file or create a new geometry. 
If you run the empty case, an error analysis will be done.
"""

from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
#from petsc4py import PETSc
import time
import os
plt.close("all")
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = 10

# define single frequency or broadband range
#Freq = np.arange(1,351,1) # frequency range, 351 max
Freq = np.array([50.0,120.0,172.0,200.0,250.0,325.0]) #104,221.29865519

# define the input files
path = "geometry/"
path_scratch = "scratch/"
name = "5Circle" #name of the .geo file : 5Circle or 5Circle_res, Empty
r = str(0.4) # define the radius of the circles.
#If you use the resonator circles, you change the inner radius, not the external one (see the .geo file).

# domain boundary (should be the same as the .geo file)
AA = 8
BB = 1


c = 343.2 # speed of sound

if name == "Empty":
    meshsize = 1/30 # use the value from gmsh
    N = np.array([4,6,8,10,15,20,30,40,60,80,100,120,140,160]) # number of elements per wavelength
    Freq = c / (meshsize*N) #associated frequency

kk = 2*np.pi*Freq/c #wavenumber
omega = 2*np.pi*Freq # pulsation
NFreq = len(Freq)
rho = 1.204 # air density
Z = rho*c # air impedance
Amp = 1 # amplitude of the wave

save_figs = 0
show_field = 1
vis_mesh = 1 # to visualize the mesh
elem_order = 2 # order of the fem
# compute pressure in dB
def PressureDB(ur,ui):
    p0 = 2e-5
    return 10*(ln((ur**2+ui**2)/(p0**2))/ln(10))
 
# Create mesh and define function space
command = 'gmsh -2 -setnumber radius '+r+' '+path+name+'.geo -o '+path_scratch+name+'.msh'

print("Running Gmsh and dolfin-convert...")
gmsh_err = os.system(command)
#gmsh_err=call(["gmsh -2 -setnumber radius "+r+" "+path + name +".geo -o "+path_scratch+name+".msh"]) #circle_centered_2D,square_2D,circle_centered_2D_small
if gmsh_err:
    print("...something bad happened when the mesh was being generated...")
else:
    dolfin_err=call(["dolfin-convert", path_scratch+name+".msh", path_scratch+name+".xml"])
if dolfin_err:
    print("...something bad happened when the mesh was being generated...")
else:
    print("...mesh generated successfully!")
# Load the xml mesh
mesh=Mesh(path_scratch+name+".xml")
N_Nodes = mesh.num_vertices()
print("...and loaded into Fenics!")
print("The mesh has",N_Nodes,"nodes")

if vis_mesh==1:
    plt.figure(num=1,figsize=(12,2))
    plot(mesh,backend="matplotlib",alpha=1)
    plt.xlabel("$x$",fontsize=font)
    plt.ylabel("$y$",fontsize=font)
    plt.xlim(-0.02,AA+0.02)
    plt.ylim(-0.02,BB+0.02)
    plt.plot([0.0,AA],[0.0,0.0],'r-',lw=2.5)
    plt.plot([0.0,AA],[BB,BB],'r-',lw=2.5)
    
    plt.plot([0.0,0.0],[0.0,BB],'b-',lw=2.5)
    plt.plot([AA,AA],[0.0,BB],'g-',lw=2.5)
    
    plt.plot(7.5,0.5,'P',markersize=8,color='black',alpha=1)
    plt.plot(7.0,0.5,'P',markersize=8,color='black',alpha=1)
    plt.plot(0.5,0.5,'P',markersize=8,color='black',alpha=1)
    plt.plot(1,0.5,'P',markersize=8,color='black',alpha=1)
    
    plt.text(-0.35,0.2,"$p_{inc}$",fontsize=font+4,color='blue')
    plt.text(3.0,1.1,"Periodic BC",fontsize=font,color='red')
    plt.text(3.0,-0.35,"Periodic BC",fontsize=font,color='red')
    plt.text(8.05,0.4,"Radiation BC",fontsize=font,color='green')
    
    plt.xticks(fontsize=font)
    plt.yticks(fontsize=font)
    if save_figs == 1:
        plt.savefig('/home/philippe/Documents/ESA/Python_SC/2D/Pics/slides/2D_MeshTrans1Dir_'+ name +'.pdf', bbox_inches='tight')

tol = 1E-8
# define the right and left borders for the BC
class Right(SubDomain):
    def inside(self, x, on_boundary):
        return (near(x[0], AA, tol)) and on_boundary 
    
class Left(SubDomain):
    def inside(self, x, on_boundary):
        return (near(x[0], 0, tol)) and on_boundary

# define boundary condations between top and bottom
class PeriodicBoundaryY(SubDomain):

    # Bottom boundary is "target domain" G
    def inside(self, x, on_boundary):
        return bool(x[1] < DOLFIN_EPS and x[1] > -DOLFIN_EPS and on_boundary)

    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        y[0] = x[0]
        y[1] = x[1] - BB

pbc = PeriodicBoundaryY()
Vr = FiniteElement("Lagrange", mesh.ufl_cell(), elem_order)
Vi = FiniteElement("Lagrange", mesh.ufl_cell(), elem_order)
Vc = Vr*Vi
V = FunctionSpace(mesh,Vc,constrained_domain=pbc)
    
mf = FacetFunction("size_t", mesh) 
mf.set_all(0) # initialize the function to zero
right = Right() # instantiate it
left = Left()
# use this left half object to set values of the mesh function to 1 in the subdomain
right.mark(mf, 1)
left.mark(mf, 2)
# define a new measure ds based on this mesh function
ds = Measure("ds")[mf] #Notation dx[meshfunction] is deprecated. Please use dx(subdomain_data=meshfunction) instead

# assuming the source to be at (0,0)
x1 = np.array([0.5,0.5])
x2 = np.array([1.0,0.5])
x3 = np.array([7,0.5]) 
x4 = np.array([7.5,0.5])

# initialize post-processing variables
Px1 = np.zeros(NFreq,dtype=complex)
Px2 = np.zeros(NFreq,dtype=complex)
Px3 = np.zeros(NFreq,dtype=complex)
Px4 = np.zeros(NFreq,dtype=complex)

A = np.zeros(NFreq,dtype=complex)
B = np.zeros(NFreq,dtype=complex)
C = np.zeros(NFreq,dtype=complex)
D = np.zeros(NFreq,dtype=complex)

error_L2r = np.zeros(NFreq)
error_L2i = np.zeros(NFreq)

# start the main frequency loop
for f in range(NFreq):
    print("Computing frequency %1.1f Hz" %(Freq[f]))
    k = 2*np.pi*Freq[f]/c
    uexr = Expression('-A*cos(K*x[0])', K = k, A = Amp, degree=4) #-
    uexi = Expression('+A*sin(K*x[0])', K = k, A = Amp, degree=4) #+
    g_i = Expression('-K', K = k, degree=4) #-
    g_r = Constant(0.0)
    
    #Define variational problem
    (u_r, u_i) = TrialFunctions(V)
    (v_r, v_i) = TestFunctions(V)

    # define variationnal formulation - term by term
    # Term 1
    t1r=(inner(grad(u_r), grad(v_r))+inner(grad(u_i), grad(v_i)))*dx
    t1i=(inner(grad(u_i), grad(v_r))-inner(grad(u_r), grad(v_i)))*dx
    # Term 2
    t2r=-k**2*(inner(u_r,v_r)+inner(u_i,v_i))*dx
    t2i=-k**2*(inner(u_i,v_r)-inner(u_r,v_i))*dx
    # Term 3
    t3r = k*(inner(u_i,v_r) - inner(u_r,v_i))*ds(1)
    t3i = -k*(inner(u_r,v_r) + inner(u_i,v_i))*ds(1) 

    ar = t1r + t2r + t3r
    ai = t1i + t2i + t3i
    a = ar+ai
        
    # assemble right hand side: int(g*v_bar domega)   (g_r+ig_i)(v_r-iv_i) = g_r*v_r + g_i*v_i + i(g_i*v_r - g_r*v_i)
    Lr = (g_r*v_r - g_i*v_i)*ds(2)
    Li = (g_i*v_r + g_r*v_i)*ds(2)
    L = Lr+Li
    
    # Compute solution
    u = Function(V)
    print("Solve linear system...")
    t = time.time()
    solve(a == L, u)
    elapsed = time.time() - t
    print("Done %1.1f sec" %elapsed)
    u_i, u_r = u.split(True)
    u_r_sc = u_r - uexr
    u_i_sc = u_i - uexi
    
    Px1[f] = u_r(x1) + 1j*u_i(x1)
    Px2[f] = u_r(x2) + 1j*u_i(x2)
    Px3[f] = u_r(x3) + 1j*u_i(x3)
    Px4[f] = u_r(x4) + 1j*u_i(x4)
    
    if name == "Empty":
        error_L2r[f] = errornorm(uexr, u_r, 'L2') # check the error for the empty case
        print('error_L2 real part =', error_L2r[f])
        error_L2i[f] = errornorm(uexi, u_i, 'L2') # check the error for the empty case
        print('error_L2 imag part =', error_L2i[f])
    
    # do the plane wave decomposition
    A[f] = 1j*(np.exp(-1j*k*x1[0])*Px2[f] - np.exp(-1j*k*x2[0])*Px1[f])/(2*np.sin(k*(x1[0]-x2[0])))   
    B[f] = 1j*(np.exp(1j*k*x2[0])*Px1[f] - np.exp(1j*k*x1[0])*Px2[f])/(2*np.sin(k*(x1[0]-x2[0])))

    C[f] = 1j*(np.exp(-1j*k*x3[0])*Px4[f] - np.exp(-1j*k*x4[0])*Px3[f])/(2*np.sin(k*(x3[0]-x4[0])))
    D[f] = 1j*(np.exp(1j*k*x4[0])*Px3[f] - np.exp(1j*k*x3[0])*Px4[f])/(2*np.sin(k*(x3[0]-x4[0])))
      
    if ((show_field == 1) and (f%5==0)):
        plt.figure(f,figsize=(5.5,2))
        #plt.subplot(211) #RdBu,RdYlBu
        p = plot(u_r,backend="matplotlib",cmap=plt.cm.get_cmap('RdBu_r', 24))#,colorbar=True,show_axis='on')#,range_min=-2., range_max=2.)
        cb = plt.colorbar(p,fraction=0.024,aspect=8,format='%.1e')
        
        cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=font-2)
        
        plt.ylabel("$y$",fontsize=font)
        plt.xlabel("$x$",fontsize=font)
        plt.xticks(fontsize=font)
        plt.yticks(fontsize=font)
        plt.title("Pressure real part [Pa], Frequency = %1.1f Hz" %(Freq[f]),fontsize=font)#%(2*np.pi*Freq[f]/c),fontsize=14)
        #cax1= plt.axes([0.92, 0.18, 0.015, 0.2])
        #,cax = cax1)#

#        plt.subplot(212)
#        p = plot(PressureDB(u_r,u_i),backend="matplotlib",cmap=plt.cm.get_cmap('magma', 24))#,range_min=70., range_max=110.)
#        #p.set_cmap("magma") #YlGnBu
#        #plt.colorbar(p,fraction=0.008)#,ticks=v)#,cax = cax2) #plt.colorbar(p, fraction=0.01, pad=0.02)
#        mini = 0
#        maxi = 100
#        #p.set_clim(mini,maxi)
#        cbar = plt.colorbar(p,fraction=0.018,aspect=8,format='%1.1f')#ticks=np.linspace(mini, maxi, num=11, endpoint=True))
#        #
#        #cbar.set_clim(vmin=mini,vmax=maxi)
#        plt.ylabel("$y$",fontsize=font)
#        plt.xlabel("$x$",fontsize=font)
#        plt.xticks(fontsize=font)
#        plt.yticks(fontsize=font)
#        plt.title("Pressure [dB]",fontsize=font)#%(2*np.pi*Freq[f]/c),fontsize=14)
        
        if save_figs == 1:
            plt.savefig('/home/philippe/Documents/ESA/Python_SC/2D/Transmission/SC_1Dir/plots/ReIm'+name+str(Freq[f])+'.pdf', bbox_inches='tight')
        #plt.close(f)

if name != "Empty":    
    T = D/B 
    R = A/B
    TL = 10*np.log10(abs(T))
    # Plot transmission loss and coefficients
    plt.figure(998)
    plt.plot(Freq,TL,"b-+")
    #print(sum(TL))
    plt.grid(True)
    plt.figure(999)
    plt.plot(Freq,abs(T),"b-+")
    plt.plot(Freq,abs(R),"r-+")
    plt.plot(Freq,abs(R)**2+abs(T)**2,"g-+")
    plt.grid(True)

#path_bis = "/home/philippe/Documents/ESA/Python_SC/2D/Transmission/SC_1Dir/data_coeffs/"
#if save_figs == 1:
#    np.savetxt(path_bis + "T_coeff_"+name+".txt",T)
#    np.savetxt(path_bis + "R_coeff"+name+".txt",R)
#    np.savetxt(path_bis + "Freq_domain"+name+".txt",Freq)
    
if name == "Empty":
    abs_err_L2 = np.sqrt(error_L2r**2+error_L2r**2)
    plt.figure(1000)
    plt.loglog(N,abs_err_L2,"r-+")
    plt.xlabel("element per wavelength",fontsize=font)
    plt.ylabel("L2 error",fontsize=font)
    plt.xticks(fontsize=font)
    plt.yticks(fontsize=font)
    #path_err= "/home/philippe/Documents/ESA/Python_SC/2D/Transmission/SC_1Dir/error_data/"
    #np.savetxt(path_err + "errL2_order4_mesh30.txt",abs_err_L2)
plt.show()