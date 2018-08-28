#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 09:51:20 2018

@author: philippe

example of usage: python SC_bempp.py 450 meshes/cylinders_coarse.msh
"""
import bempp.api
from bempp.api.operators.potential import helmholtz as helmholtz_potential
import numpy as np
import time
import matplotlib.cm as cm
from matplotlib.patches import Circle
from matplotlib import pyplot as plt
import sys
import os
from scipy.sparse.linalg import gmres
print("modules imported !")

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
                        
def run_bempp_sc(freq,meshfile):

    ################# define parameters ###################
    freq = float(freq)
    omega = 2*np.pi*freq
    c = 343.0
    k = omega/c # wavenumber
    d = np.array([0., 1., 0]) # direction of the incident plane wave
    d /= np.linalg.norm(d)
    Ampli = 1.0
    p0 = 2e-5
    z_cut = 0.0

    tol_gmres = 1e-5 
    
    vis_figs = 1 # set 1 to visualize the results with matplotlib 

    # location of the source
    xs = 0.6 
    ys = -0.5
    zs = 0.5
    source = "monopole" # monopole or plane wave
    
    ############### define input pressure ################
    def u_inc(x):
        if source == "plane_wave":
            return np.exp(1j * k * (d[0]*x[0] + d[1]*x[1] + d[2]*x[2]))
        elif source == "monopole":
            r = np.sqrt((x[0]-xs)**2+(x[1]-ys)**2+(x[2]-zs)**2)
            return Ampli*np.exp(1j * k * r)/r
        else: 
            raise ValueError('not implemented yet')
        
    def dirichlet_fun(x, n, domain_index, result):
        result[0] = np.exp(1j * k * np.dot(x, d))
    
    def neumann_fun(x, n, domain_index, result):
        result[0] = 1j * k * np.dot(n, d) * np.exp(1j * k * np.dot(x, d))
     
    def dirichlet_fun_monopole(x, n, domain_index, result): 
        r = np.sqrt((x[0]-xs)**2+(x[1]-ys)**2+(x[2]-zs)**2)
        result[0] = Ampli*np.exp(1j * k * r)/r
        
    def neumann_fun_monopole(x, n, domain_index, result):
        r = np.sqrt((x[0]-xs)**2+(x[1]-ys)**2+(x[2]-zs)**2)
        result[0]= Ampli/r*(1j*k-1/r)*np.exp(1j*k*r)* (((x[0]-xs)*n[0] + (x[1]-ys)*n[1] + (x[2]-zs)*n[2])/r)


    ################# load mesh #########################################
    grid = bempp.api.import_grid(meshfile)
    space = bempp.api.function_space(grid, "P", 1)
    print("Mesh successfully loaded !")
    print("The space has {0} dofs".format(space.global_dof_count))

    ################# define BEM formulation ###################################
    identity = bempp.api.operators.boundary.sparse.identity(
        space, space, space)
    dlp = bempp.api.operators.boundary.helmholtz.double_layer(
        space, space, space, k)
    hyp = bempp.api.operators.boundary.helmholtz.hypersingular(
        space, space, space, k, use_slp=True)
    ntd = bempp.api.operators.boundary.helmholtz.osrc_ntd(space, k)
    
    burton_miller = .5 * identity - dlp - ntd * hyp
    
    if source == "plane_wave":
        dirichlet_grid_fun = bempp.api.GridFunction(space, fun=dirichlet_fun)
        neumann_grid_fun = bempp.api.GridFunction(space, fun=neumann_fun)
        rhs_fun = dirichlet_grid_fun - ntd * neumann_grid_fun
    elif source == "monopole":
        dirichlet_grid_fun_monopole = bempp.api.GridFunction(space, fun=dirichlet_fun_monopole)
        neumann_grid_fun_monopole = bempp.api.GridFunction(space, fun=neumann_fun_monopole)
        rhs_fun = dirichlet_grid_fun_monopole - ntd * neumann_grid_fun_monopole
    else: 
        raise ValueError('not implemented yet')
     
    # bem assembling    
    print("Assembling BEM operator...")
    t = time.time()
    discrete_op = burton_miller.strong_form()
    coeffs = rhs_fun.coefficients
    elapsed = time.time() - t
    print("Bem operator assembled in %1.1f sec" %elapsed)

    # solve linear system
    t = time.time()
    counter = gmres_counter()
    x, info = gmres(discrete_op, coeffs,x0 = coeffs, maxiter = 200,restrt=100,tol = tol_gmres,callback=counter)
    elapsed = time.time() - t
    print("Gmres solving time: %1.1f sec" %elapsed)
    It = counter.niter
    Residuals = np.asarray(counter.residuals)

    total_field = bempp.api.GridFunction(discrete_op, coefficients=x)
    



##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Post processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ############ compute results on a polar grid ###########
    theta = np.linspace(0, 2 * np.pi, 600)
    R = 1.45
    xR = 0.6
    yR = 0.45
    points = np.array([xR + R*np.cos(theta), yR + R*np.sin(theta), z_cut*np.ones(len(theta))])
    slp_pot_polar = helmholtz_potential.double_layer(
        space, points, k)
    
    res_polar = u_inc(points) + slp_pot_polar.evaluate(total_field)    
    u_polar = np.abs(res_polar.flat)
    u_in = u_inc(points)
    
    if vis_figs == 1:
        plt.figure(figsize=(10, 8))
        plt.polar(theta, u_polar,"r-",linewidth=2.0)
        plt.polar(theta, np.abs(u_in),"g-.",linewidth=2.0)  

    ############ compute results on grid #################
    t = time.time()
    print("Computing the external field...")
    Nx = 100 
    Ny = 100 
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
    
    
    res = u_inc(points_grid) + slp_pot.evaluate(total_field)
    u_evaluated = res.reshape((Nx, Ny))
    elapsed = time.time() - t
    print("Time: %1.1f sec" %elapsed)
    
    ############ plot results on grid #################
    if vis_figs == 1:
        fig,ax = plt.subplots(1,figsize=(10, 8))
        ax.set_aspect('equal')
        cmap = cm.magma
        
        mini = 35
        maxi = 110
        levels = np.linspace(mini,maxi,24) #maxi-mini+1
        Z = pressure_db(u_evaluated,p0)
    
        p = ax.contourf(x.reshape((Nx, Ny)), y.reshape((Nx, Ny)), Z, levels,
                     cmap=cm.get_cmap(cmap, len(levels)))
        p2 = ax.contour(x.reshape((Nx, Ny)), y.reshape((Nx, Ny)), Z, p.levels, colors='white',linewidths=0.5,linestyles='solid',alpha=0.4)
   
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
        ax.set_title("Total Pressure map, solution in plane z=%1.1f" %z_cut)
        cbar = fig.colorbar(p)
        cbar.set_label('Pressure (dB)', rotation=90)
        plt.show()
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if __name__ == "__main__":
    import resource
    if len(sys.argv) > 2:
        print("***************************************************************")
        print("************* bempp for high frequency scattering *************")
        print("***************************************************************")

        print("Running case " + str(sys.argv[2]) )
        print("Running frequency " + str(sys.argv[1]) + " Hz")
        
        run_bempp_sc(sys.argv[1],sys.argv[2])
        print("frequency " +str(sys.argv[1]) + " Hz finished")            
        used_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        print("used memory: " +str(used_mem/1000)+ " Mb")
        print("---------------------------------------------------------------")
        print("---------------------------------------------------------------")
    else:
        raise sys.exit("usage:  python " +str(sys.argv[0])+ " <frequency> <path_to_the_mesh_file.msh>")
