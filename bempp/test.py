#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 09:51:20 2018

@author: philippe
"""

# import bempp modules
import bempp.api
from bempp.api.operators.potential import helmholtz as helmholtz_potential

# import usual python librairies
import numpy as np
import time
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Circle
import sys

# import special functions
from scipy.special import spherical_jn as jn
from scipy.special import spherical_yn as yn
from scipy.special import eval_legendre
from scipy.sparse.linalg import gmres
# define analytical solution: scattering from the unit sphere 
def analytic_sol(R,a,theta,Ampli,k,N_terms):
    result = 0
    for n in range(N_terms):
        tosum = 0
        cn = -Ampli*(2*n+1)*(-1j)**(n)*(jn(n,k*a,derivative=True)/ (jn(n,k*a,derivative=True) - 1j*yn(n,k*a,derivative=True)))    
        tosum = cn*(jn(n,k*R,derivative=False) - 1j*yn(n,k*R,derivative=False))*eval_legendre(n,np.cos(theta))
        result = result + tosum
    return result

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
            



def run_bempp_sphere(freq):
    # Parameters
    #freq = 546.2197646913849
    #freq = 1638.6592940741543
    freq = float(freq)
    omega = 2*np.pi*freq
    c = 343.2
    #rho = 1.204
    k = omega/c
    Ampli = 1
    p0 = 2e-5
    vis_mesh = 0
    vis_figs = 0
    save_data = 1
    # linear solver parameters
    tol_gmres = 1e-5
        
    # define input boundary condition: plane wave in the x-direction
    def dirichlet_fun(x, n, domain_index, result):
        result[0] = np.exp(1j * k * x[0])
        
    def neumann_fun(x, n, domain_index, result):
        result[0] = 1j * k * n[0] * np.exp(1j * k * x[0])
    
    # create the grid with approximately 10 elements per wavelength
    h = 2*np.pi/(10*k)
    grid = bempp.api.shapes.sphere(h=h)
    space = bempp.api.function_space(grid, "P", 1)
    print("Mesh successfully loaded !")
    print("The space has {0} dofs".format(space.global_dof_count))
    
    if vis_mesh == 1:
        grid.plot()
     
    # define the operators to set up the bem problem   
    identity = bempp.api.operators.boundary.sparse.identity(
        space, space, space)
    dlp = bempp.api.operators.boundary.helmholtz.double_layer(
        space, space, space, k)
    hyp = bempp.api.operators.boundary.helmholtz.hypersingular(
        space, space, space, k, use_slp=True)
    ntd = bempp.api.operators.boundary.helmholtz.osrc_ntd(space, k)
    
    # define the regularized burton miller formulation
    burton_miller = .5 * identity - dlp - ntd * hyp
    
    dirichlet_grid_fun = bempp.api.GridFunction(space, fun=dirichlet_fun)
    neumann_grid_fun = bempp.api.GridFunction(space, fun=neumann_fun)
    rhs_fun = dirichlet_grid_fun - ntd * neumann_grid_fun
    
    #A = bempp.api.as_matrix(burton_miller.weak_form())
    #print("Condition number of the bem matrix = %1.2f" %np.linalg.cond(bem_mat))
      
    t = time.time()
    discrete_op = burton_miller.strong_form()
    coeffs = rhs_fun.coefficients
    elapsed = time.time() - t
    print("Assembling bem operator: %1.1f sec" %elapsed)

    # solve the linear system 
    t = time.time()
    counter = gmres_counter()
    x, info = gmres(discrete_op, coeffs, tol = tol_gmres,callback=counter)
    elapsed = time.time() - t
    print("Gmres solving time: %1.1f sec" %elapsed)
    It = counter.niter
    Residuals = np.asarray(counter.residuals)
    total_field = bempp.api.GridFunction(discrete_op, coefficients=x)
    
    if vis_figs == 1:
        plt.figure(figsize=(9, 3))
        plt.semilogy(range(0,It), Residuals,"b-*",linewidth=2.0,label="residuals",markersize=10,markeredgecolor="k")
        plt.semilogy(range(0,It), tol_gmres*np.ones(It),"k--",linewidth=1.5,label="gmres tolerance")
        plt.legend()
        plt.grid(True,which="both",ls="--",alpha=0.5)
        plt.title("Gmres convergence")
    if save_data == 1:
        np.savetxt("Residuals_gmres.txt",Residuals)
    
    ####################################################################################
    # plot polar pressure response in the xy-plane
    theta = np.linspace(0, 2 * np.pi, 200)
    xc = 0
    yc = 0
    z_cut=0
    R = 3
    
    # define the bempp solution
    points = np.array([(R-xc)*np.cos(theta), (R-yc)*np.sin(theta), z_cut*np.ones(len(theta))])
    slp_pot_polar = helmholtz_potential.double_layer(
        space, points, k)
    res_polar = np.exp(1j * k * points[0]) + slp_pot_polar.evaluate(total_field)
    u_polar = np.abs(res_polar.flat/(np.exp(1j *k * points[0])))
    
    # define the analytical solution
    u_in = np.exp(-1j * k * points[0])
    u_polar_ex = np.abs(analytic_sol(R,1,theta,Ampli,k,100) + u_in)
    
    if vis_figs == 1:
        plt.figure(figsize=(10, 8))
        plt.polar(theta, u_polar,"r-.",linewidth=2.0,label="bempp")
        plt.polar(theta, u_polar_ex,"b:",linewidth=2.0,label="analytical")
        plt.legend()
    if save_data == 1:
        np.savetxt("polar_data"+str(freq)+"Hz.txt",u_polar)
        np.savetxt("polar_data_ex"+str(freq)+"Hz.txt",u_polar_ex)
    err = np.abs(u_polar - u_polar_ex)
    print("max error polar plot = %1.3e" %np.max(err))
    # compute the numerical solution on a grid
    t = time.time()
    Nx = 200 
    Ny = 200
    xmin, xmax, ymin, ymax = [-3, 3, -3, 3]
    plot_grid = np.mgrid[xmin:xmax:Nx * 1j, ymin:ymax:Ny * 1j]
    points_grid = np.vstack((plot_grid[0].ravel(),
                        plot_grid[1].ravel(),
                        z_cut*np.ones(plot_grid[0].size)))
    u_evaluated = np.zeros(points_grid.shape[1], dtype=np.complex128)
    u_evaluated[:] = np.nan
    
    
    x, y, z = points_grid
    slp_pot = helmholtz_potential.double_layer(
        space, points_grid, k)
    res = np.exp(1j * k * points_grid[0]) + slp_pot.evaluate(total_field)
    u_evaluated = res.reshape((Nx, Ny))
    elapsed = time.time() - t
    print("Matrix-vector product: %1.1f sec" %elapsed)
    
    # visualize results
    if vis_figs == 1:
        fig,ax = plt.subplots(1,figsize=(10, 8))
        ax.set_aspect('equal')
        cmap = cm.magma
        #p = ax.imshow(pressure_db(u_evaluated.T,p0), interpolation='bilinear',extent=[xmin, xmax, ymin, ymax],
        #cmap=cmap, vmin = 80, vmax =100,origin='lower')
        mini = 80
        maxi = 100
        levels = np.linspace(mini,maxi,maxi-mini+1)
        Z = pressure_db(u_evaluated,p0)
    
        p = ax.contourf(x.reshape((Nx, Ny)), y.reshape((Nx, Ny)), Z, levels,
                     cmap=cm.get_cmap(cmap, len(levels)))
        p2 = ax.contour(x.reshape((Nx, Ny)), y.reshape((Nx, Ny)), Z, p.levels, colors='white',linewidths=0.5,linestyles='solid',alpha=0.5) #dashed, dashdot
    
        circ = Circle((0,0),1.05,color="white")
        ax.add_patch(circ)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        cbar = plt.colorbar(p,ticks=levels)
        cbar.set_label('Pressure (dB)', rotation=90)
        ax.set_title("Scattering from the unit sphere, solution in plane z=%1.1f" %z_cut)
        plt.show()
        
    if save_data == 1:
        np.savetxt("grid_data_real"+str(freq)+"Hz.txt",np.real(u_evaluated))
        np.savetxt("grid_data_imag"+str(freq)+"Hz.txt",np.imag(u_evaluated))
        #np.savetxt("grid_points_x.txt",x)
        #np.savetxt("grid_points_y.txt",y)
        #np.savetxt("grid_points_z.txt",z)
if __name__ == "__main__":
    #freq = 163.86592940741542
    if len(sys.argv) > 1:
        run_bempp_sphere(sys.argv[1])
    else:
        raise sys.exit("usage:  python " +str(sys.argv[0])+ " <frequency>")