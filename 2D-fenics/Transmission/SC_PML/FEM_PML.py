#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solves the Helmholtz transmission problem in 2D with PML's in both x and y-directions.
Simulates the acoustic scattering of a finite sonic crystal.

Set up:
Plane wave input on the left border
PML on the top, bottom and right borders. Dirichlet BC (p=0) at the end of the PML's.

The running frequency should be chosen w.r.t the PML efficiency.
"""
from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
from toolbox import *
from subprocess import call
from petsc4py import PETSc
import time
import scipy.special as sp
from scipy.integrate import simps
plt.close("all")
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = 14

#Freq  =np.array([50,75,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450])  #np.array([172.0]) 
#Freq = np.arange(50,452,2)
Freq = np.array([172.0])
NFreq = len(Freq)

c = 343.2
kk = 2*np.pi*Freq/c
rho = 1.204
Z = rho*c
omega = 2*np.pi*Freq

path = "geometry/"
name = "RectCircle9_res_new"#RectCircle9_new,RectCircleEmpty_new, RectCircle9_res_new
PML=True
err_analysis = False # set to False when not using the empty mesh !!!
IL_analysis = True
#savepath="/home/philippe/Documents/ESA/Report/fig/2D/PML/"
save=False
showfigs = 1
grid,x_full_max,y_full_max,R,L_PML = read_gmsh_param(path,name,PML) # this reads the data from the .geo file, cf toolbox.py

x_full_max=x_full_max+L_PML
y_full_max=y_full_max+L_PML
x_full_min = 0.0 
y_full_min = -L_PML # -2
order_fem = 1


Amp = 1
Sigma_pml = 1.0*c#5 # PML parameter for n=1 linear 150Hz
#Sigma_pml = 1.0*c #1.6 bermudez modified 150Hz
# define the limits of the physical domain
x_real_min = 0.0
x_real_max = x_full_max - L_PML
y_real_max = y_full_max - L_PML
y_real_min = 0.0
    

# Create mesh and define function space
print("Running Gmsh and dolfin-convert...")
gmsh_err=call(["gmsh", "-2", path + name +".geo", "-o", path+"tmp.msh"]) #circle_centered_2D,square_2D,circle_centered_2D_small
if gmsh_err:
    print("...something bad happened when the mesh was being generated...")
else:
    dolfin_err=call(["dolfin-convert", path+"tmp.msh", path+"tmp.xml"])
if dolfin_err:
    print("...something bad happened when the mesh was being generated...")
else:
    print("...mesh generated successfully!")
# Load the xml mesh
mesh=Mesh(path+"tmp.xml")
print("...and loaded into Fenics!")

N_Nodes = mesh.num_vertices()
print("Number of nodes: ", N_Nodes)


#plot_mesh(mesh,x_full_max,y_full_max,x_full_min,y_full_min,L_PML,save,savepath,name)

tol = 1E-8
# define the limits of the FULL domain physical+PML
class Up(SubDomain):
    def inside(self, x, on_boundary):
        return (near(x[1], y_full_max, tol)) and on_boundary

class Down(SubDomain):
    def inside(self, x, on_boundary):
        return (near(x[1], y_full_min, tol)) and on_boundary

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return (near(x[0], x_full_max, tol)) and on_boundary # 10 for 6 res, 7 for 3 res, 5 for one res
    
class Left(SubDomain):
    def inside(self, x, on_boundary):
        return (near(x[0], x_full_min, tol)) and on_boundary #and (x[1]<=6.0) and (x[1]>=2.0)

Vr = FiniteElement("Lagrange", mesh.ufl_cell(), order_fem)
Vi = FiniteElement("Lagrange", mesh.ufl_cell(), order_fem)
Vc = Vr*Vi
V = FunctionSpace(mesh,Vc)

bc_R = DirichletBC(V, Constant((0.0, 0.0)), Right())
bc_D = DirichletBC(V, Constant((0.0, 0.0)), Down())
bc_U = DirichletBC(V, Constant((0.0, 0.0)), Up())  
bc = [bc_R,bc_D,bc_U]

mf = FacetFunction("size_t", mesh) 
mf.set_all(0) # initialize the function to zero
    
right = Right() # instantiate it
left = Left()
top = Up()
bottom = Down()
# use this lefthalf object to set values of the mesh function to 1 in the subdomain
right.mark(mf, 1)
left.mark(mf, 2)
top.mark(mf,3)
bottom.mark(mf,4)
## define a new measure ds based on this mesh function
ds = Measure("ds")[mf] #Notation dx[meshfunction] is deprecated. Please use dx(subdomain_data=meshfunction) instead
#ds = ds(subdomain_data = mf)

class Omega_physical(SubDomain):
    def inside(self, x, on_boundary):
        return ((x[0] <= x_real_max) and (x[1] <= y_real_max) and (x[1] >= y_real_min) and (x[0] >= x_real_min))

class Omega_Record_IL(SubDomain):
    def inside(self, x, on_boundary):
        return ((x[0] <= x_real_max-1) and (x[1] <= y_real_max-1) and (x[1] >= y_real_min+1) and (x[0] >= x_real_max-3))    

class Omega_Record_RL(SubDomain):
    def inside(self, x, on_boundary):
        return ((x[0] <= 1.4) and (x[1] <= y_real_max-1) and (x[1] >= y_real_min+1) and (x[0] >= x_real_min)) 
    
domains = CellFunction("size_t", mesh)

subdomain_physical = Omega_physical() 
subdomain_physical.mark(domains, 5)
physical_mesh = SubMesh(mesh, domains, 5)
physical_space_r = FunctionSpace(physical_mesh, Vr)
physical_space_i = FunctionSpace(physical_mesh, Vi)

subdomain_IL = Omega_Record_IL() 
subdomain_IL.mark(domains, 6)
IL_mesh = SubMesh(mesh,domains,6)
IL_space_r = FunctionSpace(IL_mesh, Vr)
IL_space_i = FunctionSpace(IL_mesh, Vi) 

subdomain_RL = Omega_Record_RL() 
subdomain_RL.mark(domains, 7)
RL_mesh = SubMesh(mesh,domains,7)
RL_space_r = FunctionSpace(RL_mesh, Vr)
RL_space_i = FunctionSpace(RL_mesh, Vi) 

domains.set_all(0)

error_L2r = np.zeros(NFreq)
error_L2i = np.zeros(NFreq)
IL_integrated = np.zeros(NFreq)
RL_integrated = np.zeros(NFreq)

x0 = np.array([8.0,5.0]) #+5 #8 for 6res, 5 for 3 res ########7
Px = np.zeros(NFreq,dtype=complex)

for f in range(NFreq):
    
    print("Computing frequency %1.1f Hz" %(Freq[f]))
    k = 2*np.pi*Freq[f]/c
    omega = k*c
    uexi = Expression('A*cos(K*x[0])', K = k, A = Amp, degree=2) #np.exp(1j * k * x[0])
    uexr = Expression('A*sin(K*x[0])', K = k, A = Amp, degree=2)
    g_i = Expression('-K', K = k, degree=2)
    g_r = Constant(0.0)
    
    #Define variational problem
    (u_r, u_i) = TrialFunctions(V)
    (v_r, v_i) = TestFunctions(V)

    ################# weak form with PML ######### #LimPhy_min, LimPhy_max, LimDom_min,LimDom_max,Sigma, omega,element,vect
    inv_alphaX = inv_alpha_Expression(LimPhy_min = x_real_min, LimPhy_max = x_real_max, LimDom_min = x_full_min, LimDom_max = x_full_max,Sigma = Sigma_pml, omega = omega,element=V.ufl_element(),vect=0)
    alphaX = alpha_Expression(LimPhy_min = x_real_min, LimPhy_max = x_real_max, LimDom_min = x_full_min, LimDom_max = x_full_max,Sigma = Sigma_pml, omega = omega,element=V.ufl_element(),vect=0) 
    # alpha[0] partie reelle, alpha[1] partie imag
    inv_alphaY = inv_alpha_Expression(LimPhy_min = y_real_min, LimPhy_max = y_real_max, LimDom_min = y_full_min, LimDom_max = y_full_max,Sigma = Sigma_pml, omega = omega,element=V.ufl_element(),vect=1)
    alphaY = alpha_Expression(LimPhy_min = y_real_min, LimPhy_max = y_real_max, LimDom_min = y_full_min, LimDom_max = y_full_max,Sigma = Sigma_pml, omega = omega,element=V.ufl_element(),vect=1)
    
    dxR = alphaX[0]*alphaY[0] + alphaX[1]*alphaY[1] / (alphaY[0]**2+alphaY[1]**2)
    dxI = alphaX[1]*alphaY[0] - alphaY[1]*alphaX[0] / (alphaY[0]**2+alphaY[1]**2)
    dyR = alphaY[0]*alphaX[0] + alphaY[1]*alphaX[1] / (alphaX[0]**2+alphaX[1]**2)
    dyI = alphaY[1]*alphaX[0] - alphaX[1]*alphaY[0] / (alphaX[0]**2+alphaX[1]**2)
    
    xyR = inv_alphaX[0]*inv_alphaY[0] - inv_alphaX[1]*inv_alphaY[1]
    xyI = inv_alphaX[1]*inv_alphaY[0] + inv_alphaX[0]*inv_alphaY[1]
    
    t1r_X, t1i_X = Assemble_PML_terms(dxR,dxI,grad(u_r)[0],grad(u_i)[0],grad(v_r)[0],grad(v_i)[0])
    t1r_Y, t1i_Y = Assemble_PML_terms(dyR,dyI,grad(u_r)[1],grad(u_i)[1],grad(v_r)[1],grad(v_i)[1])
    t2r, t2i = Assemble_PML_terms(xyR,xyI,u_r,u_i,v_r,v_i) 

    ar = t1r_X + t1r_Y -k**2*t2r
    ai = t1i_X + t1i_Y -k**2*t2i
    a = ar+ai
    ####################################################
    # assemble right hand side
    # int(g*v_bar domega)   (g_r+ig_i)(v_r-iv_i) = g_r*v_r + g_i*v_i + i(g_i*v_r - g_r*v_i)
    Lr = (g_r*v_r + g_i*v_i)*ds(2)
    Li = (g_i*v_r - g_r*v_i)*ds(2) 
    L = Lr+Li
    
    # Compute solution
    u = Function(V)
    print("Solve linear system...")
    t = time.time()
    solve(a == L, u,bc)
    elapsed = time.time() - t
    print("Done %1.1f sec" %elapsed)
    u_i, u_r = u.split(True)
    
    if IL_analysis == True:

        ur_IL = project(u_r, IL_space_r).vector()
        ui_IL = project(u_i, IL_space_i).vector()
        IL_db = 20*np.log10(1/2e-5) - 10*np.log10((ur_IL*ur_IL+ui_IL*ui_IL)/(2e-5**2))
        IL_integrated[f] = np.sum(IL_db)/len(IL_db)
        
        ur_RL = project(u_r, RL_space_r).vector()
        ui_RL = project(u_i, RL_space_i).vector()
        RL_db = 20*np.log10(1/2e-5) - 10*np.log10((ur_RL*ur_RL+ui_RL*ui_RL)/(2e-5**2))
        RL_integrated[f] = np.sum(RL_db)/len(RL_db)
        print('Integrated IL (dB)', IL_integrated[f])
        print('Integrated RL (dB)', RL_integrated[f])
        
        del ur_IL, ui_IL, ur_RL, ui_RL
    
    if err_analysis == True:
        
        u_real = project(u_r, physical_space_r)
        u_real_ex = project(uexr, physical_space_r)
        error_L2r[f] = errornorm(u_real, u_real_ex, 'L2')
        print('error_L2 real part =', error_L2r[f])
        
        u_imag = project(u_i, physical_space_i) #physical_space_i
        u_imag_ex = project(uexi, physical_space_i)
        error_L2i[f] = errornorm(u_imag,u_imag_ex, 'L2')    
        print('error_L2 imag part =', error_L2i[f])
    
    Px[f] = u_r(x0) + 1j*u_i(x0)
    
    if showfigs == 1:
#        plt.figure(2*f,figsize=(6,6)) #12,10
#        p = plot(u_r,backend="matplotlib")#,colorbar=True,show_axis='on')#,range_min=-1.1, range_max=1.1)  #u_imag,u_i  #u_imag-u_imag_ex
#        p.set_cmap("RdBu_r")
#        plt.colorbar(p)
#        plt.ylabel("$y$",fontsize=14)
#        plt.xlabel("$x$",fontsize=14)
#        plt.xlim([0,10])
#        plt.ylim([0,10])
        
#        plt.figure(f+2,figsize=(6,6)) #12,10
#        p = plot(u_imag_ex,backend="matplotlib",colorbar=True,show_axis='on')#,range_min=-1.1, range_max=1.1)
#        p.set_cmap("RdBu_r")
#        plt.colorbar(p)
#        plt.ylabel("$y$",fontsize=14)
#        plt.xlabel("$x$",fontsize=14)

        plt.figure(2*f+1,figsize=(5,5)) #12,10
        #u_imag-u_imag_ex
        mini = 35
        maxi = 110
        
        p_db = Pressure_DB(u_r,u_i)
        p = plot(p_db,backend="matplotlib",cmap=plt.cm.get_cmap('magma', 32),levels=np.linspace(mini, maxi,(maxi-mini+1),endpoint=True))#,colorbar=True,show_axis='on')#,range_min=-1.1, range_max=1.1)  #u_imag,u_i  #u_imag-u_imag_ex
        #p.set_cmap("magma")
        #plt.colorbar(p)
        cb = plt.colorbar(p,fraction=0.0425,aspect=22,format='%.1f',ticks=np.linspace(mini, maxi,(int((maxi-mini)/5)),endpoint=True))
        cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=14)
        
        plt.ylabel("$y$",fontsize=16)
        plt.xlabel("$x$",fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlim([0.0,10.0])
        plt.ylim([0.0,10.0])
        if save == True:
            plt.savefig(savepath+'PmapdB_'+ name +str(int(Freq[f]))+'.pdf', bbox_inches='tight')
            plt.savefig(savepath+'PmapdB_'+ name +str(int(Freq[f]))+'.png', bbox_inches='tight')
        
        
IL = -20*np.log10(abs(Px)/2e-5) + 20*np.log10(1/2e-5)
print("IL at x0: ",IL)
err = np.sqrt(error_L2r**2+error_L2i**2)
#print("err",err)

# save IL, RL and frequency range
#np.savetxt("/home/philippe/Documents/ESA/Python_SC/2D/Transmission/SC_PML/data/IL_"+name+".txt",IL_integrated)
#np.savetxt("/home/philippe/Documents/ESA/Python_SC/2D/Transmission/SC_PML/data/RL_"+name+".txt",RL_integrated)
#np.savetxt("/home/philippe/Documents/ESA/Python_SC_dev/Data_PML/Frequency_range.txt",Freq)

# save error data
#np.savetxt("/home/philippe/Documents/ESA/Python_SC_dev/Data_err_PML/errL2_delta4.txt",err)

#plt.figure(999)
#plt.plot(Freq,IL_integrated,"b--*",lw=1,markersize=9)
#plt.plot(Freq,RL_integrated,"r--o",lw=1,markersize=6)
#plt.ylim(0,30)

plt.show()