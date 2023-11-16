# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 22:42:07 2021

@author: MD Soyabbir Abu Hanif
"""

import math
import numpy as np
import scipy.special as sc
import matplotlib.pyplot as plt

r = 200*10**(-6)
a = 200*(10**-6) #radius of the particle
rho_o = 997 # density of host fluid
rho_i = 2100 # density of silica/scatterer
rho_cap = 2.106
vo = 1497 # outside speed (water)
vi = 5968 # inside speed (silica)
mu0_s = 30.9*(10**9)
mu_pp = 0.000891 # mu = mup - iw mu_pp
eta = 0 # shear viscosity
att_o = 0.023*(10**-12)
att_i = 2.6*(10**-22)
pi= 3.14159
f=np.logspace(-4, 5.825 , 500)
kc = 2*pi*(f/vo) + 1j*att_o*f**2
kcp = 2*pi*f/vi + 1j*att_i*f**2
ks = (1+1j)*(np.sqrt(rho_o*2*pi*f/(2*mu_pp)))
ksp =(np.sqrt(rho_i))*(2*pi*f)/(np.sqrt(mu0_s))
n = 0

yc = kc*a
ycp = kcp*a
ys = ys1=ys2=ys3=ys4 = ks*a
ysp = ksp*a

for n in range (3): #n = 0 for monopole , n = 1 for dipole , n = 2 for quadrupole
        
    hnpyc = np.sqrt(math.pi/(2*yc))* ((n/yc)*sc.hankel1(n+.5,yc)-sc.hankel1(n+1.5,yc)) #h'1(yc)
    hnys = np.sqrt(math.pi/(2*ys))*sc.hankel1(n+.5, ys)  #h1(ys)
    jnyc = np.sqrt(math.pi/(2*yc))* sc.jv(n+.5, yc) #j1(yc)
    jnys = np.sqrt(math.pi/(2*ys))* sc.jv(n+.5, ys) #jn(ys)
    jnpys = np.sqrt(math.pi/(2*ys))* ((n/ys)*sc.jv(n+.5,ys)-sc.jv(n+1.5,ys)) #jn(ys')
    jnpyc  =  np.sqrt(math.pi/(2*yc))* ((n/yc)*sc.jv(n+.5,yc)-sc.jv(n+1.5,yc)) #j'1(yc)
    jnpycp = np.sqrt(math.pi/(2*ycp))* ((n/ycp)*sc.jv(n+.5,ycp)-sc.jv(n+1.5,ycp)) #j'1(yc')
    jnysp = np.sqrt(math.pi/(2*ysp))* sc.jv(n+.5, ysp) #j1(ys')
    jnpysp =  np.sqrt(math.pi/(2*ysp))* ((n/ysp)*sc.jv(n+.5,ysp)-sc.jv(n+1.5,ysp)) #j'1(ys')
    jnycp = np.sqrt(math.pi/(2*ycp))* sc.jv(n+.5, ycp) #j'1(yc')
    hnyc = np.sqrt(math.pi/(2*ys))*sc.hankel1(n+.5, ys)  #h1(yc)
    hnpys = np.sqrt(math.pi/(2*ys))* ((n/ys)*sc.hankel1(n+.5,ys)-sc.hankel1(n+1.5,ys)) #h'1(ys)
        
    if n == 0 :
            
        a = (rho_cap/ysp**2)*((ysp**2 *jnycp) +(4*ycp*jnpycp))*yc*jnpyc
        b = (1/ys**2)*((ys**2*jnyc)+(4*yc*jnpyc))*ycp*jnpycp
        c = (1/ys**2)*((ys**2*hnyc)+(4*yc*hnpyc))*ycp*jnpycp
        d = yc*hnpyc *(rho_cap/ysp**2)*((ysp**2*jnycp)+(4*ycp*jnpycp))
        T0cc = (a-b)/(c-d)
        
    elif n ==1:
        a1 = 1j*yc**3 *((hnys)-(ys*hnpys))*(rho_cap-1)
        b1 = 3*(((4*rho_cap)-7)*hnys +((1+2*rho_cap)*ys*hnpys))
        T1cc  = -(a1/b1)
            
        a2 = (rho_cap-1)*yc
        b2 = ((4*rho_cap-7)*hnys) + ((1+2*rho_cap)*ys*hnpys)
        T1cs  = (a2/b2)
            
        a3 = 2*(rho_cap-1)*yc**2
        b3 = ys*(((4*rho_cap-7)*hnys)+((2*rho_cap+1)*ys*hnpys))
        T1sc  = -(a3/b3)
            
        a4 = ((4*rho_cap-7)*jnys) + ((2*rho_cap +1)*ys*jnpys)
        b4 = ((4*rho_cap -7)*hnys)+ ((2*rho_cap+1)*ys*hnpys)
        T1ss  = -(a4/b4)
    else :
            
        a5 = 2j*yc**5 *((ys*hnpys)-(2*hnys))
        b5 = 135*((3*hnys)+(ys*hnpys))
        T2cc = a5/b5
            
        a6 = yc**2
        b6 = 9*((3*hnys)+(ys*hnpys))
        T2cs = (a6/b6)
            
        a7 = 2*yc**2
        b7 = 3*ys *((3*hnys)+(ys*hnpys))
        T2sc = -(a7/b7)
            
        a8 = (3*jnys)+(ys*jnpys)
        b8 = (3*hnys)+(ys*hnpys)
        T2ss = a8/b8
        
        
        
        
#scattering cross section for incident compressional wave:


    

cross0 = ((1/yc)**2) * (T0cc)**2

cross1 = ((1/3*yc**2)*T1cc**2) + ( (2/(yc*ys))*(T1cs*T1cs))
cross2 = ((1/5*yc**2)*T2cc**2) + ( (6/(yc*ys))*(T2cs*T1cs))
crossSection = 4*(cross0+cross1+cross2)


plt.figure(1)

plt.plot(ys, T0cc, 'k')
plt.xlabel( ' $Y_s$' )
plt.ylabel('$T_0 ^{cc} $')
plt.title (' $ T_0^ {cc}$')
plt.grid()


plt.figure(2)
plt.plot(ys1, T1cc, 'k')
plt.xlabel( ' $Y_s$' )
plt.ylabel('$T_1 ^{cc} $')
plt.title (' $ T_1^ {cc}$')

plt.grid()


plt.figure(3)
plt.plot(ys,T1cs,'k' )
#plt.xlim()
plt.xlabel( ' $Y_s$' )
plt.ylabel('$T_1 ^{cs} $')
plt.title (' $T_1^ {cs}$')
plt.grid()

plt.figure(4)
plt.plot(ys1, T2cc, 'k')
plt.xlabel( ' $Y_s$' )
plt.ylabel('$T_2 ^{cc} $')
plt.title (' $ T_2^ {cc}$')

plt.grid()

plt.figure(4)
plt.plot(ys1, T2cs, 'k')
plt.xlabel( ' $Y_s$' )
plt.ylabel('$T_2 ^{cs} $')
plt.title (' $ T_2^ {cs}$')

plt.grid()

plt.figure(5)
plt.plot(ys1, crossSection, 'k')
plt.xlabel( ' $Y_s$' )
plt.ylabel('Cross Section ')
plt.title('Cross Section')
plt.grid()

plt.tight_layout()
plt.show()

    