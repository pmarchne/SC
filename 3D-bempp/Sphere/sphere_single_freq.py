#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computes the sphere sound-hard scattering problem fro a single frequency. Plot the polar plot pressure response.
the usage is for example: 
$ python3 sphere_single_freq.py 120"
"""
# import bempp modules
import bempp.api
from bempp.api.operators.potential import helmholtz as helmholtz_potential
# import usual python librairies
bempp.api.global_parameters.assembly.boundary_operator_assembly_type = 'hmat' #dense,hmat
bempp.api.global_parameters.assembly.potential_operator_assembly_type = 'hmat'
#bempp.api.global_parameters.assembly.enable_interpolation_for_oscillatory_kernels
bempp.api.global_parameters.hmat.admissibility = 'weak'
bempp.api.global_parameters.hmat.eps = 0.001
bempp.api.global_parameters.hmat.max_block_size = 10000

import numpy as np
import time
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Circle
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = 16
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
#    path = "/home/philippe/Documents/ESA/Python_SC/bempp/"
    freq = float(freq)
    omega = 2*np.pi*freq
    c = 343.2
    k = omega/c
    Ampli = 1
    p0 = 2e-5
    vis_mesh = 0
    # linear solver parameters
    tol_gmres = 1e-4
        
    # define input boundary condition: plane wave in the x-direction
    def dirichlet_fun(x, n, domain_index, result):
        result[0] = np.exp(1j * k * x[0])
        
    def neumann_fun(x, n, domain_index, result):
        result[0] = 1j * k * n[0] * np.exp(1j * k * x[0])
    
    # create the grid with approximately 10 elements per wavelength
    h = 2*np.pi/(10*k)
    grid = bempp.api.shapes.sphere(h=h)
    space = bempp.api.function_space(grid, "P", 1) #P,1
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

    eta = 0.03*1j/k
    burton_miller_old = .5 * identity - dlp - eta*hyp
    
    dirichlet_grid_fun = bempp.api.GridFunction(space, fun=dirichlet_fun)
    neumann_grid_fun = bempp.api.GridFunction(space, fun=neumann_fun)
    rhs_fun = dirichlet_grid_fun - ntd * neumann_grid_fun

    # compute eigenvalues and plot them
    #t = time.time()
    #full_mat = bempp.api.as_matrix(burton_miller.weak_form())
    #elapsed_hm = time.time() - t
    #full_mat_old = bempp.api.as_matrix(burton_miller_old.weak_form())
    #eig = np.linalg.eigvals(full_mat)
    #eig_old = np.linalg.eigvals(full_mat_old)
    
    # bem assembling    
    print("Assembling bem operator...")
    t = time.time()
    discrete_op = burton_miller.strong_form()
    coeffs = rhs_fun.coefficients
    elapsed_hm = time.time() - t
    print("Assembling bem operator: %1.1f sec" %elapsed_hm)

    # solve the linear system 
    t = time.time()
    counter = gmres_counter()
    x, info = gmres(discrete_op, coeffs,x0 = coeffs, maxiter = 100,tol = tol_gmres,callback=counter)
    elapsed_gmres = time.time() - t
    print("Gmres solving time: %1.1f sec" %elapsed_gmres)
    It = counter.niter
    #Residuals = np.asarray(counter.residuals)
    print("Solved in",str(It),"Iterations")
    total_field = bempp.api.GridFunction(discrete_op, coefficients=x)
        
    ####################################################################################
    # plot polar pressure response in the xy-plane
    #data_VAOne_BEM = np.loadtxt("/home/philippe/Documents/ESA/Python_SC/bempp/data_pp/Sphere/VAOne_results/Data_Sphere_FMM_546_R3.txt")
    #data_VAOne_FMM = np.loadtxt("/home/philippe/Documents/ESA/Python_SC/bempp/data_pp/Sphere/VAOne_results/Data_Sphere_FMM_546_R3.txt")
    #data_VAOne_FMM = 2e-5*10**(data_VAOne_FMM/20)
    #data_VAOne_BEM = 2e-5*10**(data_VAOne_BEM/20)

    theta = np.linspace(0, 2 * np.pi, 361)
    theta_ex = np.linspace(0, 2 * np.pi, 60)
    #theta_VAOne = np.linspace(np.pi/2, np.pi/2+2*np.pi, len(data_VAOne_BEM))
    xc = 0
    yc = 0
    z_cut=0
    R = 3
    
    # define the bempp solution
    points = np.array([(R-xc)*np.cos(theta), (R-yc)*np.sin(theta), z_cut*np.ones(len(theta))])
    points_ex = np.array([(R-xc)*np.cos(theta_ex), (R-yc)*np.sin(theta_ex), z_cut*np.ones(len(theta_ex))])
    slp_pot_polar = helmholtz_potential.double_layer(
        space, points, k)
    res_polar = np.exp(1j * k * points[0]) + slp_pot_polar.evaluate(total_field)
    u_polar = np.abs(res_polar.flat/(np.exp(1j *k * points[0])))
    
    # define the analytical solution
    #u_in = np.exp(-1j * k * points[0])
    u_in_ex = np.exp(-1j * k * points_ex[0])
    u_polar_ex = np.abs(analytic_sol(R,1,theta_ex,Ampli,k,100) + u_in_ex)
    
    fig2 = plt.figure(2,figsize=(7,7))
    ax = fig2.add_subplot(111,projection='polar')
   
    plt.plot(theta, pressure_db(u_polar,p0),"g--",linewidth=2.0,label="Bempp")
    plt.plot(theta_ex, pressure_db(u_polar_ex,p0),"ro",markersize=6,label="Analytical")
    #plt.polar(theta_VAOne, data_VAOne_BEM,"g-",linewidth=2.0,label="VAOne")
    #plt.plot(theta_VAOne, data_VAOne_FMM,"b:",linewidth=2.0,label="VAOne 546Hz")
    ax.tick_params(labelsize=font)
    ax.set_rmin(84)
    ax.set_rmax(98)
    ax.legend(fontsize=font-2)
    ax.grid(True)
    #plt.savefig("/home/philippe/Documents/ESA/Python_SC/bempp/data_pp/Sphere/plots/polar.pdf", bbox_inches='tight')

    plt.show()
    
if __name__ == "__main__":
    #freq = 163.86592940741542
    import resource
    if len(sys.argv) > 1:
        print("***************************************************************")
        print("************* bempp for high frequency scattering *************")
        print("***************************************************************")
        tmp = float(sys.argv[1])
        print("Running frequency " +str(tmp) + " Hz")
        run_bempp_sphere(tmp)
        print("frequency " +str(tmp) + " Hz finished")            
        used_mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        print("used memory: " +str(used_mem/1000)+ " Mb")
        print("---------------------------------------------------------------")
        print("---------------------------------------------------------------")
    else:
        raise sys.exit("usage:  python " +str(sys.argv[0])+ " <frequency>")
