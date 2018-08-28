#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
- Plot the direct and reciprocal lattice for square and hexagonal 2d periodic media. You can change the variable "periodicity" from 'HEX' to 'SQ'.
- Plot the 1D configuration
"""
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import rc
#import matplotlib.colors as colors
plt.close("all")
font = 10
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# select the type of periodicity: SQ or HEX
periodicity = 'HEX'
# The two real space lattice vectors
a1=np.array([1.0,0.0])
if (periodicity == 'SQ'):
    phi = np.pi/2.0 # np.pi/2 for sq periodicity
elif (periodicity == 'HEX'):
    phi = np.pi/3.0
else: 
    raise ValueError('not implemented yet')

a2 = np.array([np.cos(phi),np.sin(phi)])
# Vectors orthogonal to a1 and a2 (for plotting)
a1n=np.array([a1[1],-a1[0]])/np.sqrt(np.dot(a1,a1))
a2n=np.array([a2[1],-a2[0]])/np.sqrt(np.dot(a2,a2))

# The two reciprocal space lattice vectors
b1=2.0*np.pi*np.array([a2[1],-a2[0]])/(a1[0]*a2[1]-a1[1]*a2[0])
b2=2.0*np.pi*np.array([-a1[1],a1[0]])/(a1[0]*a2[1]-a1[1]*a2[0])

# 2N x 2N lattice points
N=10
nv=np.arange(-N,N)
mv=np.arange(-N,N)
# x-co-ordinates of the lattice points
Xv=np.array([[i*a1[0]+j*a2[0] for i in nv] for j in mv])
Yv=np.array([[i*a1[1]+j*a2[1] for i in nv] for j in mv])
# For a lattice of points [-Nmax,Nmax]x[-Nmax,Nmax]


# Plot range
xmax=0.5*N*min(np.sqrt(np.dot(a1,a1)),np.sqrt(np.dot(a2,a2)))

# Points of reciprocal space lattice
nv=np.arange(-N,N+1)
mv=np.arange(-N,N+1)
Kxv=np.array([[i*b1[0]+j*b2[0] for i in nv] for j in mv])
Kyv=np.array([[i*b1[1]+j*b2[1] for i in nv] for j in mv])

# Equations for the lines enclosing the unit cell

# Lines through the mid-point of a1
Xp1=[0.5*a1[0]+0.2*i*xmax*a2[0] for i in np.arange(-10,10)]
Yp1=[0.5*a1[1]+0.2*i*xmax*a2[1] for i in np.arange(-10,10)]
Xp2=[-0.5*a1[0]+0.2*i*xmax*a2[0] for i in np.arange(-10,10)]
Yp2=[-0.5*a1[1]+0.2*i*xmax*a2[1] for i in np.arange(-10,10)]

# Lines through the mid-point of a2
Xp3=[0.5*a2[0]+0.2*i*xmax*a1[0] for i in np.arange(-10,10)]
Yp3=[0.5*a2[1]+0.2*i*xmax*a1[1] for i in np.arange(-10,10)]
Xp4=[-0.5*a2[0]+0.2*i*xmax*a1[0] for i in np.arange(-10,10)]
Yp4=[-0.5*a2[1]+0.2*i*xmax*a1[1] for i in np.arange(-10,10)]

# Vectors to nearest neighbours
v1=b1
v2=b2
v3=-b1
v4=-b2
v5=b2+b1
v6=-b1-b2

# Vectors orthogonal to v1,v2,v3,v4,v5 and v6
v1n=np.array([v1[1],-v1[0]])
v2n=np.array([v2[1],-v2[0]])
v3n=np.array([v3[1],-v3[0]])
v4n=np.array([v4[1],-v4[0]])
v5n=np.array([v5[1],-v5[0]])
v6n=np.array([v6[1],-v6[0]])

# Lines along v1,v2,v3,v4,v5 and v6
l1=[0.5*v1+i*v1n for i in np.linspace(-10.0,10.0,10)]
l2=[0.5*v2+i*v2n for i in np.linspace(-10.0,10.0,10)]
l3=[0.5*v3+i*v3n for i in np.linspace(-10.0,10.0,10)]
l4=[0.5*v4+i*v4n for i in np.linspace(-10.0,10.0,10)]
l5=[0.5*v5+i*v5n for i in np.linspace(-10.0,10.0,10)]
l6=[0.5*v6+i*v6n for i in np.linspace(-10.0,10.0,10)]

############ Brillouin zone ###########################
# Consider a closed loop in the first Brillouin zone consisting of three points
P1=[0,0]
P2=0.5*(b1+b2)
P3=0.5*b2[1]*(b1+2*b2)/(b1[1]+2*b2[1])
# Draw the closed path
p1x=P1[0]+(P2[0]-P1[0])*np.linspace(0.0,1.0,5)
p1y=P1[1]+(P2[1]-P1[1])*np.linspace(0.0,1.0,5)
p2x=P2[0]+(P3[0]-P2[0])*np.linspace(0.0,1.0,5)
p2y=P2[1]+(P3[1]-P2[1])*np.linspace(0.0,1.0,5)
p3x=P3[0]+(P1[0]-P3[0])*np.linspace(0.0,1.0,5)
p3y=P3[1]+(P1[1]-P3[1])*np.linspace(0.0,1.0,5)
########################################################


xd=max(a1[0],a2[0])
yd=max(a1[1],a2[1])

kxd=max(b1[0],b2[0])
kyd=max(b1[1],b2[1])

xlim = 1.5
ylim = xlim
plt.figure(num=1, figsize=(5,2.5))
plt.subplot(121)
plt.plot(Xv.flatten(),Yv.flatten(),'o',markersize=24,color='black',alpha=0.2)
plt.xlim(-0.2*N*xd,0.2*N*xd)
plt.ylim(-0.2*N*yd,0.2*N*yd)
plt.xlabel("$x$",fontsize=font)
plt.ylabel("$y$",x=0.7,fontsize=font)
plt.xticks(fontsize=font)
plt.yticks(fontsize=font)
plt.title("Direct space",fontsize=font+1)
#plt.arrow(0,0,a1[0],a1[1],head_width=.025,head_length=.2,length_includes_head=True,color='red')
#plt.arrow(0,0,a2[0],a2[1],head_width=.025,head_length=.2,length_includes_head=True,color='red')
plt.text(1.1*a1[0],1.1*a1[1],"${\\bf a}_{1}$",color='red',fontsize=font)
plt.text(1.1*a2[0],1.1*a2[1],"${\\bf a}_{2}$",color='red',fontsize=font)


plt.plot(Xp1,Yp1,'--',color='black',alpha=1,linewidth=1.5)
#plt.plot(Xp1+np.linalg.norm(a1),Yp1,'k--')
plt.plot(Xp2,Yp2,'--',color='black',alpha=1,linewidth=1.5)
#plt.plot(Xp2-np.linalg.norm(a1),Yp2,'k--')
plt.plot(Xp3,Yp3,'--',color='black',alpha=1,linewidth=1.5)
#plt.plot(Xp3,Yp3+np.linalg.norm(a2[1]),'k--')
plt.plot(Xp4,Yp4,'--',color='black',alpha=1,linewidth=1.5)
#plt.plot(Xp4,Yp4-np.linalg.norm(a2[1]),'k--')
plt.arrow(0,0,a1[0],a1[1],color='r',head_width=0.1,length_includes_head=True)
plt.arrow(0,0,a2[0],a2[1],color='r',head_width=0.1,length_includes_head=True)
plt.xlim(-xlim,xlim)
plt.ylim(-ylim,ylim)
#plt.axes().set_aspect('equal')
plt.tight_layout()

plt.subplot(122)
plt.plot(Kxv.flatten(),Kyv.flatten(),'o',markersize=24,color='black',alpha=0.2)
plt.xlim(-0.2*N*kxd,0.2*N*kxd)
plt.ylim(-0.2*N*kyd,0.2*N*kyd)
plt.xlabel("$K_{x}$",fontsize=font)
plt.ylabel("$K_{y}$",x=0.7,fontsize=font)
plt.xticks(fontsize=font)
plt.yticks(fontsize=font)
plt.title("Reciprocal space",fontsize=font+1)
plt.arrow(0,0,b1[0],b1[1],color='blue',head_width=0.1*2*np.pi,length_includes_head=True)
plt.arrow(0,0,b2[0],b2[1],color='blue',head_width=0.1*2*np.pi,length_includes_head=True)
plt.text(1.0*b1[0],0.9*b1[1],"${\\bf b}_{1}$",color='blue',fontsize=font)
plt.text(1.0*b2[0]+0.5,0.95*b2[1],"${\\bf b}_{2}$",color='blue',fontsize=font)
#plt.axes().set_aspect('equal')

plt.plot([i[0] for i in l1],[i[1] for i in l1],'--',color='black',alpha=1,linewidth=1.5)
plt.plot([i[0] for i in l2],[i[1] for i in l2],'--',color='black',alpha=1,linewidth=1.5)
plt.plot([i[0] for i in l3],[i[1] for i in l3],'--',color='black',alpha=1,linewidth=1.5)
plt.plot([i[0] for i in l4],[i[1] for i in l4],'--',color='black',alpha=1,linewidth=1.5)


plt.text(P1[0]-1.1,P1[1]-1.0,"$\Gamma$",fontsize=font)

if (periodicity == 'SQ'):
    plt.plot(p1x,0*p1x,'-',color='deepskyblue')
    plt.plot(0*p1x+np.pi,p1x,'-',color='deepskyblue')
    plt.plot(p3y,p3y,'-',color='deepskyblue')

    plt.text(P2[0]+0.2,-1.35,"$X$",fontsize=font)
    plt.text(P2[0]+0.2,P2[0]-1.35,"$M$",fontsize=font)
elif (periodicity == 'HEX'):
    plt.plot([i[0] for i in l5],[i[1] for i in l5],'--',color='black',alpha=1,linewidth=1.5)
    plt.plot([i[0] for i in l6],[i[1] for i in l6],'--',color='black',alpha=1,linewidth=1.5)
    plt.plot(p1x,p1y,'-',color='deepskyblue')
    plt.plot(p2x,p2y,'-',color='deepskyblue')
    plt.plot(p3x,p3y,'-',color='deepskyblue')
 
    plt.text(P2[0],P2[1]+0.5,"$M$",fontsize=font)
    plt.text(P3[0],P3[1]+0.5,"$K$",fontsize=font)
else: 
    print('not implemented yet') 

###################################################
plt.xlim(-xlim*2*np.pi,xlim*2*np.pi)
plt.ylim(-ylim*2*np.pi,ylim*2*np.pi)
plt.xticks([-3*np.pi,-2*np.pi,-np.pi,0,np.pi,2*np.pi,3*np.pi],["$-3\pi$","$-2\pi$","$-\pi$","$0$","$\pi$","$2\pi$","$3\pi$"],fontsize=font)
plt.yticks([-3*np.pi,-2*np.pi,-np.pi,0,np.pi,2*np.pi,3*np.pi],["$-3\pi$","$-2\pi$","$-\pi$","$0$","$\pi$","$2\pi$","$3\pi$"],fontsize=font)
plt.tight_layout()

#if (periodicity == 'SQ'):
    #plt.savefig('/home/philippe/Documents/ESA/Python_SC/2D/IBZ_sq2.pdf', bbox_inches='tight')
#elif (periodicity == 'HEX'):
    #plt.savefig('/home/philippe/Documents/ESA/Python_SC/2D/IBZ_hex.pdf', bbox_inches='tight')
#else: 
    #print('not implemented yet') 

# infinite medium
width = 0.5
#plt.figure(num=2, figsize=(5,2.5))
fig2, axs = plt.subplots(1, figsize=(7/3.0,1.5))
font = 11
medium1 = np.array([2,2,2,2,2,2,2])
medium2 = np.array([2,2,2,2,2,2,2])
#ind = np.array([-1.5,-1,-0.5,0,0.5,1,1.5])
ind = np.array([-3,-2,-1,0,1,2,3])
# We can set the number of bins with the `bins` kwarg
d1 = 0.3
d2 = 0.7
axs.bar(ind - width,medium1,d1, color='k',alpha=0.3)#,width)
axs.bar(ind ,medium2, d2, color='b',alpha=0.4)

xtick = np.array([-d2/2-d1,-d2/2,d2/2])
axs.set_xticks(xtick)# + width/2)
#axs.yaxis.set_tick_params(labelsize=font)
axs.set_xticklabels(('$-h_1$','$0$','$h_2$'),fontsize=font)
axs.set_ylim(0,1)
axs.set_xlim(-d2/2-d1,d2/2)
#axs.set_title('Unit cell',fontsize=font+2)
axs.text(-0.6, 0.6, r'$\rho_1, c_1$')
axs.text(-0.1, 0.6, r'$\rho_2, c_2$')

axs.spines['left'].set_position('center')
axs.spines['bottom'].set_position('center')

# Eliminate upper and right axes
axs.spines['right'].set_color('none')
axs.spines['left'].set_visible(False)
axs.spines['top'].set_color('none')

# Show ticks in the left and lower axes only
axs.xaxis.set_ticks_position('bottom')

axs.yaxis.set_major_locator(plt.NullLocator())
axs.get_yaxis().set_ticklabels([])
#plt.savefig('/home/philippe/Documents/ESA/Python_SC/2D/Infinite_1D.pdf', bbox_inches='tight')
# finite medium: 3 unit cells
width = 0.5
fig3, ax3 = plt.subplots(1, figsize=(7,1.5))
font = 11
medium1 = np.array([2,2,2,2])
medium2 = np.array([2,2,2,2])
#ind = np.array([-1.5,-1,-0.5,0,0.5,1,1.5])
ind = np.array([-1,0,1,2])
# We can set the number of bins with the `bins` kwarg
d1 = 0.3
d2 = 0.7
ax3.bar(ind - width,medium1,d1, color='k',alpha=0.3)#,width)
ax3.bar(ind ,medium2, d2, color='b',alpha=0.4)

xtick = np.array([-3*d2/2-d1,-d2/2-d1,-d2/2,d2/2,d2/2+d1,3*d2/2+d1])
ax3.set_xticks(xtick)# + width/2)
#axs.yaxis.set_tick_params(labelsize=font)
ax3.set_xticklabels(('$-h_2-h_1$','$-h_1$','0','$h_2$','$h_2+h_1$','$2h_2$'),fontsize=font)
ax3.set_ylim(0,1)
ax3.set_xlim(-1.5*d1-1.5*d2-0.01,1.5*d1+1.5*d2+0.01)
#axs.set_title('Unit cell',fontsize=font+2)

ax3.spines['left'].set_position('center')
ax3.spines['bottom'].set_position('center')

# Eliminate upper and right axes
ax3.spines['right'].set_color('none')
ax3.spines['top'].set_color('none')

# Show ticks in the left and lower axes only
ax3.xaxis.set_ticks_position('bottom')
ax3.yaxis.set_major_locator(plt.NullLocator())
ax3.get_yaxis().set_ticklabels([])

ax3.vlines(x=0.5, ymin=0, ymax=2, linewidth=1.5,linestyle='--',color='k')
ax3.vlines(x=-0.5, ymin=0, ymax=2, linewidth=1.5,linestyle='--',color='k')
#ax3.vlines(x=1.5, ymin=0, ymax=2, linewidth=2.5,linestyle='-',color='r')
#ax3.vlines(x=-1.5, ymin=0, ymax=2, linewidth=2.5,linestyle='-',color='r')

plt.show()
#plt.savefig('/home/philippe/Documents/ESA/Python_SC/2D/finite3_1D.pdf', bbox_inches='tight')