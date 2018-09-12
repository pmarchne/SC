#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute the Band gap structure for a 2d arbitrary geometry unit-cell.
It solves the Helmoltz equation with Bloch periodic boundary conditions.

Iso-frequency countours, band gap structure along the IBZ and Bloch eigenmodes can be plotted. Turn the associated variables to 1 depending on what you want.

This routine is for geometries with square periodicity
For the mesh file: If you want to change the geometry, Try to use the same format as the other mesh files e.g Geom/SQ/Circle.geoCircle.geo
You might need to change the gridsize as well in the corresponding .geo file.
If you change the lattice constant on the .geo file, you need to do change it in the code as well.

For the resonator case, please note that you need to extend the IBZ to fully characterize the medium.

For the iso-frequency countours, the mesh file "IsoFreq_Empty.geo" defines the First Brillouin Zone.
"""

from dolfin import *
import numpy as np
import numpy.linalg as la

import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = 14

from subprocess import call
import tqdm as tqdm
#import sys
#plt.close("all")
path="Geom/SQ/"
path_scratch = "Geom/scratch/"
name = "Circle" # name of the .geo file e.g Circle, Circle_res, SQ or Circle_vaoneComp
#path_save = "/home/philippe/Documents/ESA/Report/fig/2D/BG/SQ/"
path_save = "/home/philippe/Documents/ESA/Report/fig/2D/BG/SQ/Res/"

N_eigen = 24 # number of eigenvalues to solve
pnts = 40 # 1/3 of the number of points along the IBZ

View_mesh = 1 # to view the mesh
plot_isofreq = 0 # to plot the isofrequency contours 
plot_eigenmodes = 1 # to plot the eigenmodes
plot_BG = 1 # to plot the band gap structure
save_figs = 0 # turn to 1 to save the figures, you need to specify your path_save

if plot_isofreq == 1:
    plot_eigenmodes = 0 
    plot_BG = 0

# Density and speed of sound in air - host medium
rhoh=1.204
ch=343.2

a = 1 # lattice constant
# The lattice vectors
a1=a*np.array([0,1])
a2=a*np.array([1,0])
# Reciprocal space lattice
b1=(2.0*np.pi)*np.array([a2[1],-a2[0]])/(a1[0]*a2[1]-a2[0]*a1[1])
b2=(2.0*np.pi)*np.array([a1[1],-a1[0]])/(a1[1]*a2[0]-a2[1]*a1[0])

if sum(a1*b2)==0 and sum(a2*b1)==0 and sum(a1*b1)==2*np.pi and sum(a2*b2)==2*np.pi:
    print("Reciprocal basis sucessfully generated")
else:
    print("Error in the direct or reciprocal basis")
    
# Set up mesh for 2D problem
# vertices of unit cell (in anticlockwise order), first vertex is at the bottom right corner.
# The mesh is constructed so that the periodic boundary conditions are between the side V[0]->V[1] and V[2]->V[3]
V=[[],[],[],[]]
V[0]=0.5*(-a1+a2) 
#V[0] = np.array([ 1.0, 0.0])
V[1]=0.5*(a1+a2)
#V[1] = np.array([ 1.0, 1.0])
V[2]=0.5*(a1-a2)
#V[2] = np.array([ 0.0, 1.0])
V[3]=-0.5*(a1+a2)
#V[3] = np.array([ 0.0, 0.0])

# Run gmsh and dolfin-convert to get the xml mesh file
print("Running Gmsh and dolfin-convert...")

gmsh_err=call(["gmsh", "-2", path+name+".geo", "-o", path_scratch+name+".msh"])
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
print("...and loaded into Fenics!")

if View_mesh == 1:
    plt.figure(num=0,figsize=(3,3))
    plot(mesh,backend="matplotlib")
    #plt.title("Mesh of unit cell",fontsize=font)
    plt.xlabel("$x/a$",fontsize=font+2)
    plt.ylabel("$y/a$",fontsize=font+2)
    step = 0.01
    plt.xlim(0.5*(-a1[1]+a2[1])-step,0.5*(a1[1]-a2[1])+step)
    #plt.xlim(0-step,1+step)
    plt.ylim(0.5*(-a1[0]-a2[0])-step,0.5*(a1[0]+a2[0])+step)
    #plt.ylim(0-step,1+step)
    plt.plot([V[0][0],V[1][0]],[V[0][1],V[1][1]],'b-',lw=2)
    plt.plot([V[1][0],V[2][0]],[V[1][1],V[2][1]],'r-',lw=2)
    plt.plot([V[3][0],V[0][0]],[V[3][1],V[0][1]],'r-',lw=2)
    plt.plot([V[2][0],V[3][0]],[V[2][1],V[3][1]],'b-',lw=2)
    plt.xticks(fontsize=font)
    plt.yticks(fontsize=font)
    if save_figs == 1:
        plt.savefig(path_save+'2D_Mesh'+ name +'.pdf', bbox_inches='tight')
# define the boundaries of the domain
tol = 1e-5
def Left_boundary(x):
    return near(x[0], -a1[1]/2.0, tol)
def Right_boundary(x):
    return near(x[0], a1[1]/2.0, tol)
def Up_boundary(x):   
    return near(x[1], a2[0]/2.0, tol)
def Down_boundary(x):
    return near(x[1], -a2[0]/2.0, tol)

# define the periodic boundary conditions
class PeriodicBoundary(SubDomain):
    
    def __init__(self, kx = 1., ky = 1.):
        SubDomain.__init__(self)
        self.kx = kx
        self.ky = ky
    
    # Left boundary is the "target domain"
    def inside(self, x, on_boundary):
        # return True if on left or bottom boundary AND NOT on one of the two corners (0, 1) and (1, 0)
        return bool((Left_boundary(x) or Down_boundary(x)) and 
                (not ((near(x[0], 0.5*(a1[0]+a2[0])) and near(x[1], -0.5*(a1[1]+a2[1]))) or 
                        (near(x[0], -0.5*(a1[0]+a2[0])) and near(x[1], 0.5*(a1[1]+a2[1])))) 
                  )and on_boundary)

    def map(self, x, y):
        # The top-most boundary point is identified with the lower-most one
        if near(x[0], 0.5*(a1[0]-a2[0])) and near(x[1], -0.5*(a1[1]-a2[1])):
            y[0] = x[0] - (a1[0]-a2[0])
            y[1] = x[1] + (a1[1]-a2[1])
        # The upper right boundary is shifted by a_1
        elif Up_boundary(x):
            y[0] = x[0] #- a1[0]
            y[1] = x[1] - a1[1]
        # Otherwise shift by a_2
        else:   # near(x[1], 1)
            y[0] = x[0] - a2[0]
            y[1] = x[1] #+ a2[1]


# Solver: creates and solves the variational formulation
def esol(Kx,Ky):
    K=Constant((Kx,Ky))
    # complex fields, double the function space
    Vr = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    Vi = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    Vc = Vr*Vi
    Vz = FunctionSpace(mesh,Vc,constrained_domain=PeriodicBoundary(kx = Kx, ky = Ky))
    # define real and imaginary parts for the test and trial function
    (u_r, u_i) = TrialFunctions(Vz)
    (v_r, v_i) = TestFunctions(Vz)
    # Define weak form - term by term
    # Term 1
    t1r=inner(grad(u_r), grad(v_r))-inner(grad(u_i), grad(v_i))
    t1i=inner(grad(u_r), grad(v_i))+inner(grad(u_i), grad(v_r))
    # Term 2
    t2r=v_i*inner(K,grad(u_r)) + v_r*inner(K,grad(u_i)) - u_i*inner(K,grad(v_r)) - u_r*inner(K,grad(v_i))
    t2i=-v_r*inner(K,grad(u_r)) + v_i*inner(K,grad(u_i)) + u_r*inner(K,grad(v_r)) - u_i*inner(K,grad(v_i)) 
    # Term 3
    t3r=inner(K,K)*(u_r*v_r-u_i*v_i) 
    t3i=inner(K,K)*(u_r*v_i+u_i*v_r) 

    # Sum terms
    ar = t1r+t2r+t3r
    ai = t1i+t2i+t3i
    L = Constant(0.0)*(v_r+v_i)*dx # not useful, just there for the assembly
    
    # Overlap integral between test and trial functions
    m = (u_r*v_r - u_i*v_i + u_r*v_i + u_i*v_r)*dx

    A,_ = assemble_system((ar+ai)*dx, L)#, bc)
    B = assemble(m)
    B = B/(ch**2)
    eigensolver = SLEPcEigenSolver(as_backend_type(A), as_backend_type(B))
    prm = eigensolver.parameters
    info(prm, True)
    eigensolver.parameters['spectrum'] = 'smallest real' #other options: 'largest imaginary', 'largest real','smallest magnitude'

    eigensolver.solve(N_eigen)

    return eigensolver


if plot_isofreq == 1:
    # definition of all the points from the brillouin zone for the isofrequency contours
    mesh_file_plot = path+"IsoFreq_Empty" # load the brillouin zone mesh
    gmsh_err=call(["gmsh", "-2", mesh_file_plot+".geo", "-o", path_scratch+"tmp.msh"])
    dolfin_err2=call(["dolfin-convert", path_scratch+"tmp.msh", path_scratch+"tmp.xml"])
    mesh_BZ=Mesh(path_scratch+"tmp.xml")
    xy=mesh_BZ.coordinates()
    triangles= mesh_BZ.cells()
    Kx = xy[:,0]
    Ky = xy[:,1]
    K = list(np.transpose(np.array([Kx,Ky])))
    sol=[esol(i[0],i[1]) for i in tqdm.tqdm(K)] # solve the eigenvalue problem for all the wavenumbers in the BZ
    
else:
    # definition of the path around the IBZ #
    Ks1=[np.array([np.pi*i/a,0.]) for i in np.linspace(0.001,1.0,pnts)]
    Ks1bis=[np.array([0.,np.pi*i/a]) for i in np.linspace(0.001,1.0,pnts)] # for non-symmetric geometries...
    Ks2=[np.array([np.pi/a,np.pi*i/a]) for i in np.linspace(0.001,1.0,pnts)]
    Ks2bis=[np.array([np.pi*i/a,np.pi/a]) for i in np.linspace(0.001,1.0,pnts)]
    Ks3=[np.array([np.pi*i/a,np.pi*i/a]) for i in np.linspace(1.0,0.001,pnts)]
    K=Ks1+Ks2+Ks3 
    #K=Ks1bis+Ks2bis+Ks3 # to be used for the resonant case and non-symmetric geometries
    sol=[esol(i[0],i[1]) for i in tqdm.tqdm(K)] # solve the eigenvalue problem for all the wavenumbers in the IBZ

################## plot the results #############################
if plot_isofreq == 1:
    for j in range(12):
        plt.figure(num=j+2,figsize=(4,4))
        z = [np.sqrt(np.abs(i.get_eigenpair(2*j)[0]))/(2*np.pi) for i in sol]
        Nlines  = 13
        levels = np.linspace(np.min(z),np.max(z),Nlines)
        CS = plt.tricontour(Kx,Ky,triangles,z,colors='k',levels=levels,linewidths=1,linestyles="--")
        im = plt.tricontourf(Kx,Ky,triangles,z, cmap='plasma',alpha=1,levels=levels) #,extent=[-np.pi, np.pi, -np.pi, np.pi]
        plt.clabel(CS, inline=1, fontsize=font-2,fmt='%1.1f')
        cbar = plt.colorbar(im,fraction=0.046, pad=0.08)
        plt.xticks([-np.pi,-np.pi/(2*a),0,np.pi/(2*a),np.pi],["$-\pi/a$","$-\pi/2a$","$0$","$\pi/2a$","$\pi/a$"],fontsize=font)
        plt.ylabel("$K_y$",fontsize=font)
        plt.xlabel("$K_x$",fontsize=font)
        plt.yticks([-np.pi/a,-np.pi/(2*a),0,np.pi/(2*a),np.pi/a],["$-\pi/a$","$-\pi/2a$","$0$","$\pi/2a$","$\pi/a$"],fontsize=font)
        cbar.ax.set_ylabel('Frequency [Hz]',rotation = 90,fontsize=font)
        cbar.ax.tick_params(labelsize=font)
        plt.axis('image')
        if j==0:
            color_IBZ = "grey"
            plt.annotate("$\Gamma$", xy=(0, 0),xytext=(-0.4, -0.1),fontsize=font,color=color_IBZ) #(-0.2, 0)
            plt.annotate("X", xy=(np.pi/a, 0),xytext=(np.pi/a+0.06, 0),fontsize=font,color=color_IBZ)
            plt.annotate("M", xy=(np.pi/a, np.pi/a),xytext=(np.pi/a+0.06, np.pi/a-0.1),fontsize=font,color=color_IBZ)
            plt.plot([0, np.pi/a], [0, 0], color=color_IBZ, linestyle='-', linewidth=1.5)
            plt.plot([np.pi/a, np.pi/a], [0, np.pi/a], color=color_IBZ, linestyle='-', linewidth=1.5)
            plt.plot([np.pi/a, 0], [np.pi/a, 0], color=color_IBZ, linestyle='-', linewidth=1.5)
            
            # to be used for the resonant case and non-symmetric geometries 
            plt.annotate("X'", xy=(0,np.pi/a),xytext=(-0.5,np.pi/a-0.5),fontsize=font,color=color_IBZ)
            plt.plot([0, 0], [0, np.pi/a], color=color_IBZ, linestyle='--', linewidth=1.5)
            plt.plot([0, np.pi/a], [np.pi/a-0.02, np.pi/a-0.02], color=color_IBZ, linestyle='--', linewidth=1.5)
            
        #plt.title("Equi-frequency contours - Mode " + str(j),fontsize=font)
        if save_figs ==1:
            plt.savefig(path_save+name+'2D_KmapEig_'+ str(j) +'.pdf', bbox_inches='tight')
        
# initialization of arrays to plot the band gaps 
N_gaps = 14 # estimate an upper limit of the number of BG       
BG_min = np.zeros(N_gaps) 
BG_max = np.zeros(N_gaps)
GamX = pnts
full = 3*pnts
k=0
if plot_BG == 1:
    plt.figure(num=1,figsize=(4,6))
   
    for j in range(0,N_eigen-4,2):
        plt.plot([np.sqrt(i.get_eigenpair(j)[0])/(2*np.pi) for i in sol],linewidth=2, linestyle='--') # *0.5*ch*1E-6/(s*np.pi) # /(2*np.pi)
        BG_min[k] = np.max([np.sqrt(i.get_eigenpair(j)[0])/(2*np.pi) for i in sol[0:full]])
        BG_max[k] = np.min([np.sqrt(i.get_eigenpair(j+2)[0])/(2*np.pi) for i in sol[0:full]])
        k = k+1
        if (    np.max([np.sqrt(i.get_eigenpair(j)[0])/(2*np.pi) for i in sol[0:full]]) - np.min([np.sqrt(i.get_eigenpair(j+2)[0])/(2*np.pi) for i in sol[0:full]])   <=0):
            plt.axhspan(np.max([np.sqrt(i.get_eigenpair(j)[0])/(2*np.pi) for i in sol[0:full]]),np.min([np.sqrt(i.get_eigenpair(j+2)[0])/(2*np.pi) for i in sol[0:full]]),alpha=0.3, color = 'grey')
    plt.grid(True)
    plt.xticks([0,pnts,2*pnts,3*pnts],["$\Gamma$","$X$","$M$","$\Gamma$"],fontsize=font)
    plt.ylabel("Frequency [Hz]",fontsize=font+2)
    plt.xlabel("Bloch Wavevector",fontsize=font+2)
    plt.yticks(fontsize=font)
    plt.ylim(0,500)
    plt.xlim(0,3*pnts)
    
    if save_figs ==1:
        plt.savefig(path_save+'2D_BG'+ name +'.pdf', bbox_inches='tight')
##########################################################################################################
    
if plot_eigenmodes == 1:
    n_plot = 3# number of eigenmodes to plot
    # Plot pressure field
    for n in range(n_plot): # loop over the numbers of eigenmodes to plots
        N=0# Point on K-axis e.g 29
        Vr = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
        Vi = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
        Vc = Vr*Vi
        Vz = FunctionSpace(mesh,Vc,constrained_domain=PeriodicBoundary())
        pressure = Function(Vz)
        eig=pressure.vector()
        eig[:]=sol[N].get_eigenpair(2*n)[2] #solution au point (Kx,Ky)
    
        plt.figure(num=n+3,figsize=(3,3))
        #plt.title("Pressure Eigenmode %i: Kx={1:.1f}, Ky={2:.1f}".format(n,K[N][0],K[N][1]) %n,fontsize=font)
        p = plot(pressure[0],backend="matplotlib",cmap=plt.cm.get_cmap('RdBu_r', 32))#,colorbar=True,show_axis='on')#, vmin = -1, vmax = 1)
        #p.set_cmap("plasma") #magma

        plt.xlabel("$x/a$",fontsize=font+2)
        plt.ylabel("$y/a$",fontsize=font+2)
        plt.xticks(fontsize=font)
        plt.yticks(fontsize=font)
        cbar = plt.colorbar(p, fraction=0.046, pad=0.04)
        NN = len(cbar.ax.get_yticklabels())
        ticks = np.around(np.linspace(-1,1,NN),decimals=2)
        cbar.ax.set_yticklabels(ticks, fontsize=font, visible=False)  # vertically oriented colorbar
        if save_figs == 1:
            plt.savefig(path_save+'2D_Eig'+str(n)+'_'+ name +'.pdf', bbox_inches='tight')

plt.show()
