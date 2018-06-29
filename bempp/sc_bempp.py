#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 09:51:20 2018

@author: philippe
"""
#import pyamg
import bempp.api
from bempp.api.operators.potential import helmholtz as helmholtz_potential
import numpy as np
import time
#import matplotlib.cm as cm
#from matplotlib.patches import Circle
#from matplotlib import pyplot as plt
import sys
from scipy.sparse.linalg import gmres
#plt.close("all")

def pressure_db(u,p0):
    return 10*np.log10(np.abs(u)**2/p0**2)

class gmres_counter(object):
    def __init__(self, disp=True):
        self._disp = disp
        self.niter = 0
        self.residuals = []
    def __call__(self, rk=None):
        self.niter += 1
        self.residuals.append(rk)
        if self._disp:
            print('iteration %3i: residual = %s' % (self.niter, str(rk)))
            
#class fgmres_counter(object):
#    def __init__(self, disp=True):
#        self._disp = disp
#        self.niter = 0
#    def __call__(self,xk=None):
#        self.niter += 1
#        if self._disp:
#            print('iteration %i' % (self.niter))
                        
def run_bempp_sc(freq):
    #freq = 430
    freq = float(freq)
    omega = 2*np.pi*freq
    c = 343.0
    #rho = 1.25
    k = omega/c
    d = np.array([0., 1., 0]) # direction of the incident plane wave
    d /= np.linalg.norm(d)
    #direction = 1
    Ampli = 1.0
    p0 = 2e-5
    vis_mesh = 0
    z_cut = 0.0#0.0 # longueur du cylindre / 2 , mettre a zeroooo
    
    tol_gmres = 1e-5 # tolerance for the linear solver
    
    vis_figs = 0
    save_data = 0
    
    xs = 0.6#0.6
    ys = -0.5#-0.5
    zs = 0.5
    source = "monopole"
    
    def u_inc(x):
        if source == "plane_wave":
            return np.exp(1j * k * (d[0]*x[0] + d[1]*x[1] + d[2]*x[2]))
        elif source == "monopole":
            r = np.sqrt((x[0]-xs)**2+(x[1]-ys)**2+(x[2]-zs)**2)
            return Ampli*np.exp(1j * k * r)/r#/(4*np.pi*r)
        else: 
            raise ValueError('not implemented yet')
        
    def dirichlet_fun(x, n, domain_index, result):
        result[0] = np.exp(1j * k * np.dot(x, d))
    
    def neumann_fun(x, n, domain_index, result):
        result[0] = 1j * k * np.dot(n, d) * np.exp(1j * k * np.dot(x, d))#np.exp(1j * k * x[0]), n[0]
     
    def dirichlet_fun_monopole(x, n, domain_index, result): 
        r = np.sqrt((x[0]-xs)**2+(x[1]-ys)**2+(x[2]-zs)**2)
        result[0] = Ampli*np.exp(1j * k * r)/r
        
    def neumann_fun_monopole(x, n, domain_index, result):
        r = np.sqrt((x[0]-xs)**2+(x[1]-ys)**2+(x[2]-zs)**2)
        result[0]= Ampli/r*(1j*k-1/r)*np.exp(1j*k*r)* (((x[0]-xs)*n[0] + (x[1]-ys)*n[1] + (x[2]-zs)*n[2])/r)
    
    path = "/home/philippe/Documents/ESA/Python_SC/bempp/meshes/"
    name = "cylinders3D_coarse" # use the mesh "cylinders3D_coarse" to test the code
    grid = bempp.api.import_grid(path+name+ ".msh")
    space = bempp.api.function_space(grid, "P", 1)
    print("Mesh successfully loaded !")
    print("The space has {0} dofs".format(space.global_dof_count))
    
    if vis_mesh == 1:
        grid.plot()
    
    identity = bempp.api.operators.boundary.sparse.identity(
        space, space, space)
    dlp = bempp.api.operators.boundary.helmholtz.double_layer(
        space, space, space, k)
    hyp = bempp.api.operators.boundary.helmholtz.hypersingular(
        space, space, space, k, use_slp=True)
    ntd = bempp.api.operators.boundary.helmholtz.osrc_ntd(space, k)
    
    burton_miller = .5 * identity - dlp - ntd * hyp
    
    
    #dirichlet_grid_fun.plot()
    #dirichlet_grid_fun_monopole.plot()
    #neumann_grid_fun_monopole.plot()
    #neumann_grid_fun.plot()
    dirichlet_grid_fun = bempp.api.GridFunction(space, fun=dirichlet_fun)
    neumann_grid_fun = bempp.api.GridFunction(space, fun=neumann_fun)
    if source == "plane_wave":
        rhs_fun = dirichlet_grid_fun - ntd * neumann_grid_fun
    elif source == "monopole":
        dirichlet_grid_fun_monopole = bempp.api.GridFunction(space, fun=dirichlet_fun_monopole)
        neumann_grid_fun_monopole = bempp.api.GridFunction(space, fun=neumann_fun_monopole)
        rhs_fun = dirichlet_grid_fun_monopole - ntd * neumann_grid_fun_monopole
    else: 
        raise ValueError('not implemented yet')
     
    # bem assembling    
    print("Assembling bem operator...")
    t = time.time()
    discrete_op = burton_miller.strong_form()
    coeffs = rhs_fun.coefficients
    elapsed = time.time() - t
    print("Bem operator assembled in %1.1f sec" %elapsed)

    # solve linear system
    t = time.time()
    counter = gmres_counter()
    x, info = gmres(discrete_op, coeffs,x0 = coeffs, maxiter = 100,restrt=50,tol = tol_gmres,callback=counter)
    elapsed = time.time() - t
    print("Gmres solving time: %1.1f sec" %elapsed)
    It = counter.niter
    Residuals = np.asarray(counter.residuals)
    
#    rk=[]
#    t = time.time()
#    counter = fgmres_counter()
#    x,info = pyamg.krylov.gmres(discrete_op, coeffs, x0 = coeffs, tol = tol_gmres, restrt=None, maxiter=100,callback=counter,residuals=rk)
#    elapsed = time.time() - t
#    print("fGmres solving time: %1.1f sec" %elapsed)
#    It = counter.niter
#    Residuals_pyamggmres = np.asarray(rk)
    total_field = bempp.api.GridFunction(discrete_op, coefficients=x)
    
    if vis_figs == 1:
        plt.figure(figsize=(6, 8))
        plt.semilogy(range(0,It), Residuals,"r-*",linewidth=2.0,label="residuals",markersize=10,markeredgecolor="k")
        plt.semilogy(range(0,It), tol_gmres*np.ones(It),"k--",linewidth=2.0,label="gmres tolerance")
        plt.legend()
        plt.grid(True,which="both",ls="--",alpha=0.5)
        plt.title("Gmres convergence")
    if save_data == 1:
        np.savetxt("ResidualsSC_gmres.txt",Residuals)
    ####################################################################################
    #polar plots
    theta = np.linspace(0, 2 * np.pi, 600)
    R = 1.45
    xR = 0.6
    yR = 0.45
    #zR = 0.0
    points = np.array([xR + R*np.cos(theta), yR + R*np.sin(theta), z_cut*np.ones(len(theta))])
    slp_pot_polar = helmholtz_potential.double_layer(
        space, points, k)
    
    #r = np.sqrt((points[0]-xs)**2+(points[1]-ys)**2+(points[2]-zs)**2)
    res_polar = u_inc(points) + slp_pot_polar.evaluate(total_field)    
    u_polar = np.abs(res_polar.flat)#/u_inc(points))
    u_in = u_inc(points)
    
    if vis_figs == 1:
        plt.figure(figsize=(10, 8))
        plt.polar(theta, u_polar,"r-",linewidth=2.0)
        plt.polar(theta, np.abs(u_in),"g-.",linewidth=2.0)  
    if save_data == 1:
        np.savetxt("polar_dataSC"+str(freq)+"Hz.txt",u_polar)
    # compute insertion loss
    coor_IL = np.array([[0.6],[1.9],[0.0]])
    slp_pot_IL = helmholtz_potential.double_layer(
        space, coor_IL, k)
    res_IL = u_inc(coor_IL) + slp_pot_IL.evaluate(total_field)
    IL = float(pressure_db(u_inc(coor_IL),p0)-pressure_db(res_IL,p0))
    print("Point coordinates xyz: %1.2f %1.2f %1.2f" %(coor_IL[0] ,coor_IL[1] ,coor_IL[2]))
    print("Corresponding Insertion loss: %1.2f dB" %IL)
    ############ grid plot #################
    t = time.time()
    print("Computing matrix-vector product...")
    Nx = 300 #300
    Ny = 300 #300
    xmin, xmax, ymin, ymax = [-0.5, 1.7,-0.5,2.2]# limits of the pressure map #0.5
    plot_grid = np.mgrid[xmin:xmax:Nx*1j , ymin:ymax:Ny*1j] #* 1j]
    points_grid = np.vstack((plot_grid[0].ravel(),
                        plot_grid[1].ravel(),
                        z_cut*np.ones(plot_grid[0].size)))
    u_evaluated = np.zeros(points_grid.shape[1], dtype=np.complex128)
    u_evaluated[:] = np.nan
    
    
    x, y, z = points_grid
    slp_pot = helmholtz_potential.double_layer(
        space, points_grid, k)
    
    
    res = u_inc(points_grid) + slp_pot.evaluate(total_field) #np.exp(1j * k * points[0])
    u_evaluated = res.reshape((Nx, Ny))
    elapsed = time.time() - t
    print("Matrix-vector product: %1.1f sec" %elapsed)
    
    if vis_figs == 1:
        fig,ax = plt.subplots(1,figsize=(10, 8))
        ax.set_aspect('equal')
        cmap = cm.magma
        
        mini = 35
        maxi = 110
        levels = np.linspace(mini,maxi,maxi-mini+1)
        Z = pressure_db(u_evaluated,p0)
    
        p = ax.contourf(x.reshape((Nx, Ny)), y.reshape((Nx, Ny)), Z, levels,
                     cmap=cm.get_cmap(cmap, len(levels)))
        p2 = ax.contour(x.reshape((Nx, Ny)), y.reshape((Nx, Ny)), Z, p.levels, colors='white',linewidths=0.5,linestyles='solid',alpha=0.4) #dashed, dashdot
    
        #p = ax.imshow(pressure_db(u_evaluated.T,p0), interpolation='bilinear',extent=[xmin, xmax, ymin, ymax],cmap=cm.magma,origin='lower')#, vmin = 87, vmax =98)
   
        a = 0.3 # lattice constant
        x_coor = np.array([0.0,a,2*a,3*a,4*a])
        x_coor = np.repeat(x_coor,4)
        y_coor = np.array([0.0,0.3,0.6,0.9])
        y_coor = np.tile(y_coor,5)
        R_circ = 0.1
        
        for i in range(len(x_coor)):
            ax.add_patch(Circle((x_coor[i],y_coor[i]),R_circ*(1+0.01),color="white"))
    
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title("Pressure map, solution in plane z=%1.1f" %z_cut)
        cbar = fig.colorbar(p)
        cbar.set_label('Pressure (dB)', rotation=90)
        plt.show()
    if save_data == 1:
        np.savetxt("grid_data_realSC"+str(freq)+"Hz.txt",np.real(u_evaluated))
        np.savetxt("grid_data_imagSC"+str(freq)+"Hz.txt",np.imag(u_evaluated))
if __name__ == "__main__":
    #freq = 687.5 or 430 Hz
    if len(sys.argv) > 1:
        run_bempp_sc(sys.argv[1])
    else:
        raise sys.exit("usage:  python " +str(sys.argv[0])+ " <frequency>")
