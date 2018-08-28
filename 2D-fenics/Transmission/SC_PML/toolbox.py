#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
additional functions for the acoustic transmission problem with PML
"""
import re
from dolfin import *
import matplotlib.pyplot as plt
import numpy as np


def Pressure_DB(ur,ui):
    return (10/ln(10))*ln((ur**2+ui**2)/(2e-5**2))
        
        
def read_gmsh_param(path,name,PML):

    file = open(path+name+".geo", "r") 
    maxlines = 15
    head = [next(file) for x in range(maxlines)]

    head = ' '.join(head).split()
    del head[0]
    del head[0]

    gridsize = re.findall("\d+", head[0])
    gridsize= float(gridsize[0])/float(gridsize[1])

    a = re.findall("\d+", head[1])
    a = float(a[0])
    
    b = re.findall("\d+", head[2])
    b = float(b[0])

    r = re.findall("\d+\.\d+", head[3])
    r = float(r[0])
    
    if PML == True:
        L_PML = re.findall("\d+", head[4])
        L_PML = float(L_PML[0])
    else:
        L_PML=0.0

    print("gridsize: %1.4f" %(gridsize),"\nlength x-dir:",a,"\nlength y-dir:",b,"\ncircle_radius:",r,"\nwidth PML:",L_PML)
    return gridsize,a,b,r,L_PML


def Assemble_PML_terms(alphar,alphai,pr,pi,qr,qi):  
    Real_part = alphar*inner(pr,qr) + alphar*inner(pi,qi) - alphai*inner(pi,qr) + alphai*inner(pr,qi) 
    Imag_part = alphar*inner(pi,qr) - alphar*inner(pr,qi) + alphai*inner(pr,qr) + alphai*inner(pi,qi) 
    return Real_part*dx, Imag_part*dx

class alpha_Expression(Expression):

  def __init__(self, LimPhy_min, LimPhy_max, LimDom_min,LimDom_max,Sigma, omega,element,vect): #vect: 0 in x direction, 1 in y direction
        #SubDomain.__init__(self)
        self.LimPhy_max = LimPhy_max
        self.LimPhy_min = LimPhy_min
        self.LimDom_max = LimDom_max
        self.LimDom_min = LimDom_min
        self.Sigma = Sigma
        self.omega = omega
        self._ufl_element = element
        self.vect = vect
    
  def eval(self, value, x):
    n = 1
    delta = self.LimDom_max - self.LimPhy_max
    step = 0.01
    if (x[self.vect] <= self.LimPhy_max) and (x[self.vect] >= self.LimPhy_min):
        value[0] = 1.0
        value[1] = 0.0
    else:  
        if(x[self.vect] < self.LimPhy_min):
            #SigmaF = self.Sigma*((x[self.vect]-self.LimPhy_min)/(self.LimDom_min-self.LimPhy_min))**n
            SigmaF = self.Sigma/(delta-(-x[self.vect]-self.LimPhy_min)+step) #- self.Sigma/(delta+step) # Bermudez normalized
        else:
            #SigmaF = self.Sigma*((x[self.vect]-self.LimPhy_max)/(self.LimDom_max-self.LimPhy_max))**n 
            SigmaF = self.Sigma/(delta-(x[self.vect]-self.LimPhy_max)+step) #- self.Sigma/(delta+step) # Bermudez normalized
        value[0] = 1.0/(1+(SigmaF**2/self.omega**2)) #1.0/(1+(self.Sigma**2/k**2))
        value[1] = -(SigmaF/self.omega) / (1+(SigmaF**2/self.omega**2))
  def value_shape(self):
    return (2,)

class inv_alpha_Expression(Expression):

  def __init__(self, LimPhy_min, LimPhy_max, LimDom_min,LimDom_max,Sigma, omega,element,vect):
        #SubDomain.__init__(self)
        self.LimPhy_max = LimPhy_max
        self.LimPhy_min = LimPhy_min
        self.LimDom_max = LimDom_max
        self.LimDom_min = LimDom_min
        self.Sigma = Sigma
        self.omega = omega
        self._ufl_element = element
        self.vect = vect
    
  def eval(self, value, x):
    n=1
    delta = self.LimDom_max - self.LimPhy_max
    step = 0.01
    if (x[self.vect] <= self.LimPhy_max) and (x[self.vect] >= self.LimPhy_min):
        value[0] = 1.0
        value[1] = 0.0
    else:
        if(x[self.vect] < self.LimPhy_min):
            #SigmaF = self.Sigma*((x[self.vect]-self.LimPhy_min)/(self.LimDom_min-self.LimPhy_min))**n
            SigmaF = self.Sigma/(delta-(-x[self.vect]-self.LimPhy_min)+step) #- self.Sigma/(delta+step) # Bermudez normalized
        else:
            #SigmaF = self.Sigma*((x[self.vect]-self.LimPhy_max)/(self.LimDom_max-self.LimPhy_max))**n    
            SigmaF = self.Sigma/(delta-(x[self.vect]-self.LimPhy_max)+step) #- self.Sigma/(delta+step) # Bermudez normalized
        value[0] = 1.0
        value[1] = SigmaF/self.omega
  def value_shape(self):
    return (2,)


def plot_mesh(mesh,x_max,y_max,x_min,y_min,L,save,savepath,name):
    plt.figure(num=0,figsize=(4,5))
    plot(mesh,backend="matplotlib",alpha=0.35)
    #plt.title("Mesh of unit cell",fontsize=16)
    #plt.title("${\\rm Mesh\;of\;unit\;cell\;(indicating\;identified\;sides)}$")
    plt.xlabel("$x$",fontsize=16)
    plt.ylabel("$y$",fontsize=16)
    plt.xlim(x_min-0.1,x_max+0.1)
    plt.ylim(y_min-0.05,y_max+0.05)
    plt.plot([0.0,x_max],[y_max-L,y_max-L],'g-',lw=2.0)
    plt.plot([0.0,x_max],[y_min,y_min],'m-',lw=2.0)
    plt.plot([0.0,x_max],[y_min+L,y_min+L],'g-',lw=2.0)
    plt.plot([0.0,x_max],[y_max,y_max],'m-',lw=2.0)
    
    plt.plot([x_min,x_min],[y_min,y_max],'b-',lw=2.0)
    plt.plot([x_max-L,x_max-L],[y_min,y_max],'g-',lw=2.0)
    plt.plot([x_max,x_max],[y_min,y_max],'m-',lw=2.0)
      
    x_real_max = x_max - L
    y_real_max = y_max - L
    y_real_min = 0
    plt.plot([x_real_max-1,x_real_max-1],[y_real_min+1,y_real_max-1],'k--',lw=2.0)
    plt.plot([x_real_max-3,x_real_max-3],[y_real_min+1,y_real_max-1],'k--',lw=2.0)
    plt.plot([x_real_max-3,x_real_max-1],[y_real_min+1,y_real_min+1],'k--',lw=2.0)
    plt.plot([x_real_max-3,x_real_max-1],[y_real_max-1,y_real_max-1],'k--',lw=2.0)

    
    plt.text(-3.0,6.5,"$p_{inc}$",fontsize=18,color='blue')
    plt.text(5.0,y_max+0.2,"Dirichlet BC",fontsize=14,color='m')
    plt.text(5.0,y_max-L+0.2,"PML",fontsize=16,color='green')
    plt.text(5.0,-1.4,"PML",fontsize=16,color='green')
    plt.text(10.5,5.0,"PML",fontsize=16,color='green')
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    if save == True:
        plt.savefig(savepath+'2D_MeshTrans2Dir_'+ name +'.pdf', bbox_inches='tight')



