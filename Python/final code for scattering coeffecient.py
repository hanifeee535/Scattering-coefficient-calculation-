"""
Created on Wed Feb 17 14:14:52 2021
Title:  A python code for calculating scattering coefficient both analytical and matrix 
form and comparing them with their visual representations using
Numpy , Scipy and Matplotlib

@author: Md SOyabbir Abu Hanif 
dept. of EEE, Sylhet Engineering College, Sylhet , Bangladesh.
"""


import math
from decimal import Decimal
from decimal import *
import numpy as np

import scipy.special as sc
import matplotlib.pyplot as plt
crossC  = []
ysys  = []
t0ccarray = []
t1ccarray = []
t1csarray = []
t1scarray = []
t1ssarray = []
t2ccarray = []
t2csarray = []
t2scarray = []
t2ssarray = []


t0ccreal =[]
t1ccreal =[]
t1csreal =[]
t2ccreal =[]
t2csreal =[]

t0ccimag =[]
t1ccimag =[]
t1csimag =[]
t2ccimag =[]
t2csimag =[]


#physical properties 

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


#frequency span 
ff=np.logspace(-4,5.825, 500)
nf = len (ff)
kc = 2*pi*(ff/vo) + 1j*att_o*ff**2
kcp = 2*pi*ff/vi + 1j*att_i*ff**2
ks = (1+1j)*(np.sqrt(rho_o*2*pi*ff/(2*mu_pp)))
ksp =(np.sqrt(rho_i))*(2*pi*ff)/(np.sqrt(mu0_s))
ab  = ks*a
yc = kc*a
ys = ks*a
ycp = kcp*a
ysp = ksp*a

order_N = 2

#analytical calculation : 

def analytical ():
    TA = []
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
    ys = ks*a
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
            
            
    TA = [T0cc, T1cc , T1cs , T1sc ,T1ss,T2cc , T2cs, T2sc , T2ss]
    return TA

coeff = analytical()

t0cc = coeff[0]
t1cc = coeff[1]
t1cs = coeff[2]
t1sc = coeff[3]
t1ss = coeff[4]
t2cc = coeff[5]
t2cs = coeff[6]
t2sc = coeff[7]
t2ss = coeff[8]

        
    




#Matrix coefficent function  

def coefficient ( f):
    T = []
    
    

    a = (200*(10**-6))
    rho_i = 2100
    rho_o = 997
    rho_cap = 2.106
    vo = 1497 # outside speed (water)
    vi = 5968 # inside speed (silica)
    mu0_s = 30.9*(10**9)
    mu_pp = 0.000891 # mu = mup - iw mu_pp
    eta = 0 # shear viscosity
    att_o = 0.023*(10**-12)
    att_i = 2.6*(10**-22)
    pi= 3.14159
    kc = 2*pi*(f/vo) + 1j*att_o*f**2
    kcp = 2*pi*f/vi + 1j*att_i*f**2
    ks = (1+1j)*(np.sqrt(rho_o*2*pi*f/(2*mu_pp)))
    ksp =(np.sqrt(rho_i))*(2*pi*f)/(np.sqrt(mu0_s))
    
    
    
    ys = ks*a
    yc = kc*a
    ysp = ksp*a
    ycp = kcp*a
    
    
    for order in range (2): #m = 0 for compressional incident , m = 1 for shear incident 
        
        
        for n in range (3): #n =0 for monopole , n = 1 for dipole , no = 2 for quadropole
            
            #calculating bessel and hunnkel function and there derivetives:
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
            
            #calculating matrix components :
            
            m11 = yc*hnpyc
     
        
            m12 = n*(n+1)*hnys
    

        
            m13 = -(ycp * jnpycp)

        
            m14 = -(n*(n+1)*jnysp)

        
            m21 = hnyc
   
       
            m22 = (ys*hnpys) + hnys 

        
            m23 = -(jnycp)

        
            m24  =  -(ysp*jnpysp) + jnysp  
        
            m31 = ((ys**2- 2*n*(n+1))*hnyc + 4*yc*hnpyc)*(1/ys**2)

        
            m32 = (hnys - ys*hnpys)*2*n*(n+1)*(1/ys**2)

        
            m33  = -((ysp**2- (2*n*(n+1)))*jnycp +4*ycp*jnpycp)*(rho_cap/ysp**2)

       

            m34  = -( jnysp -ysp*jnpysp) *n*(n+1)*(2*rho_cap/ysp**2)

        
            m41 = (hnyc - (yc*hnpyc)) * (2/ys**2) 

        
            m42 = (2*ys*hnpys + (ys**2 -(2*n*(n+1))+2)*hnys)*(1/ys**2)

        
            m43 = -2* (jnycp - (ycp*jnpycp))* (rho_cap/ysp**2)
        
            m44 = - ( 2*ysp*jnpysp + ( ysp**2-(2*n*(n+1))+2)*jnysp ) * (rho_cap/ysp**2)

        
            c1 = -(yc*jnpyc)


        
            c2 = -jnyc

       
            c3 = - ( (ys**2-(2*n*(n+1))) *jnyc + 4*yc*jnpyc )* (1/ys**2)

        
            c4 =  -(jnyc - (yc*jnpyc)) *(2/ys**2)
        
            c5 = -(n*(n+1))*jnys
        
            c6 = -(ys*jnpys+jnys)
        
            c7 = -((2/ys**2) * (n*n+1)*((ys*jnpys)-jnys))
        
            c8 = -(((2*n*(n+1))-2-ys**2)*jnys -2*ys*jnpys) *(1/ys**2)
        
            M = np.array([
                [m11,m12,m13,m14],
                [m21,m22,m23,m24],
                [m31,m32,m33,m34],
                [m41,m42,m43,m44]  

                ])

            c = np.array([[c1,c5],[c2,c6],[c3,c7],[c4,c8]]) 
        
        
            col1 =  (abs(np.mean(M[:,0])))
            col2 =  (abs(np.mean(M[:,1])))
            col3 =  abs(np.mean(M[:,2]))
            col4 =  abs(np.mean(M[:,3]))
            col5 =  abs(np.mean(c[:,0]))
            col6 =  abs(np.mean(c[:,1]))
        
        
            OrderMagnitude1 = ((math.log(col1,10)))
            OrderMagnitude2 = (math.log(col2))
            OrderMagnitude3 = (math.log(col3, 10))
            OrderMagnitude4 = (math.log(col4, 10))
            OrderMagnitude5 = (math.log(col5, 10))
            OrderMagnitude6 = (math.log(col6, 10))


            M[:,0] = M[:,0] * 10**(-OrderMagnitude1)
            M[:,1] = M[:,1] * 10**(-OrderMagnitude2)
            M[:,2] = M[:,2] * 10**(-OrderMagnitude3)
            M[:,3] = M[:,3] * 10**(-OrderMagnitude4)
            c[:,0] = c[:,0] * 10**(-OrderMagnitude5)
            c[:,1] = c[:,1] * 10**(-OrderMagnitude5)
            
            if n == 0:

                M0 = np.array([
                    [m11*10**(-OrderMagnitude1) , m13*10**(-OrderMagnitude3)],
                    [m31*10**(-OrderMagnitude1) , m33*10**(-OrderMagnitude3)]
                ])



                c0 = np.array ([
                    [c1*10**(-OrderMagnitude5)],
                    [c3*10**(-OrderMagnitude5)]
                ])



                inverse=np.linalg.inv(M0)
                X0 = np.dot(M0,c0)
                



                if order == 0 :
                    T0cc = X0[0]*10**(OrderMagnitude5-OrderMagnitude1)
                    
            elif n ==1:
                inverse1 = np.linalg.inv(M)
                
                if order == 0:
                    X1 = np.dot(inverse1, c[:,0] )
                    T1cc = X1[0]*10**(OrderMagnitude5-OrderMagnitude1)
                    T1cs = X1[1]*10**(OrderMagnitude5-OrderMagnitude2)
                   
                    
                if order  == 1:
                    X2 = np.dot(inverse1, c[:,1] )
                    T1sc = X2[0]*10**(OrderMagnitude6-OrderMagnitude1)
                    #T1ss = complex((X2[1]*10**(OrderMagnitude6-OrderMagnitude2)))
                    T1ss = 0
            else:
                inverse2 = np.linalg.inv(M)
                
                if order == 0:
                    X3 = np.dot(inverse1, c[:,0] )
                    T2cc = X3[0]*10**(OrderMagnitude5-OrderMagnitude1)
                    T2cs = X3[1]*10**(OrderMagnitude5-OrderMagnitude2)
                   
                    
                if order  == 1:
                    X4 = np.dot(inverse1, c[:,1] )
                    T2sc = X4[0]*10**(OrderMagnitude6-OrderMagnitude1)
                    T2ss =0
                    #T2ss = X4[1]*10**(OrderMagnitude6-OrderMagnitude2)
                   
                   
                    

                   
            

    T = [T0cc , T1cc, T1cs , T1sc, T1ss, T2cc , T2cs , T2sc, T2ss]
    return T
            
        
               
        

        
        
            
    

#matrix equation calculation 
for kf in range (nf):

    frequency = ff[kf]

    TT = coefficient ( frequency)

    aa = TT[0]
    bb = TT[1]
    cc = TT[2]
    dd = TT[3]
    ee = TT[4]
    jj = TT[5]
    gg = TT[6]
    hh =TT[7]
    ii = TT[8]

    t0ccarray.append (aa)
    #t0ccarray[kf] = aa
    
    t1ccarray.append (bb)
    #t1ccarray[kf] = bb
    t1csarray.append (cc)
    #t1csarray[kf] = cc
    t1scarray.append (dd)
    t1ssarray.append (ee)
    t2ccarray.append (jj)
    t2csarray.append (gg)
    t2scarray.append (hh)
    t2ssarray.append (ii)
    
    t0ccreal.append(aa.real)
    t1ccreal.append(bb.real)
    t1csreal.append(cc.real)
    t2ccreal.append(jj.real)
    t2csreal.append(gg.real)
    
    t0ccimag.append(aa.imag)
    t1ccimag.append(bb.imag)
    t1csimag.append(cc.imag)
    t2ccimag.append(jj.imag)
    t2csimag.append(gg.imag)
    

#scattering cross section for incident compressional wave:



    
for y in range (nf) :
    co1 = t0ccarray[y]
    co2 = t1ccarray[y]
    co3 = t1csarray[y]
    co4 = t2ccarray[y]
    co5 = t2csarray[y]
    
    
    
    
    cross0 = (co1**2)/(yc[y])**2  

    cross1 = ((co2**2)/3*(yc[y])**2) + ( (2/(yc[y]*ys[y]))*(co3**2))
    cross2 = ((co4**2)/5*(yc[y])**2) + ( (6/(yc[y]*ys[y]))*(co5**2))
    crossSection = 4*(cross0+cross1+cross2)
    crossC.append(crossSection)
    
    
    

#separating real and imaginary part of the coefficients
    
    

    

#Plotting the graphs




plt.figure(1)

plt.plot(ab,t0ccarray,'k',  )
#plt.plot (ab , t0cc , 'b--', label = 'Analytical')
plt.xlabel( '$Y_s$' )
plt.ylabel(' $T_0 ^{cc} $')
plt.title (' $T_0^ {cc}$')
plt.grid()
plt.legend()

plt.figure(2)

plt.plot(ab,t1ccarray, 'k' ,)
#plt.plot (ab , t1cc , 'b--', label = 'Analytical')
plt.xlabel( '$Y_s$' )
plt.ylabel(' $T_1 ^{cc} $')
plt.title (' $T_1^ {cc}$')
plt.grid()
plt.legend()

plt.figure(3)

plt.plot(ab,t1csarray, 'k', )
#plt.plot (ab , t1cs , 'b--', label = 'Analytical')
plt.xlabel( '$Y_s$' )
plt.ylabel(' $T_1 ^{cs} $')
plt.title (' $T_1^ {cs}$')
plt.grid()
plt.legend()

plt.figure(4)

plt.plot(ab,t1scarray, 'b', )
#plt.plot (ab , t1sc , 'b--', label = 'Analytical')
plt.xlabel( '$Y_s$' )
plt.ylabel(' $T_1 ^{sc} $')
plt.title ('Dipole coefficient $T_1^ {sc}$')
plt.grid()
plt.legend()

plt.figure(5)

plt.plot(ab,t1ssarray, 'b', )
#plt.plot (ab , t1ss , 'b--', label = 'Analytical')
plt.xlabel( '$Y_s$' )
plt.ylabel(' $T_1 ^{ss} $')
plt.title ('Dipole coefficient $T_1^ {ss}$')
plt.grid()
plt.legend()

plt.figure(6)

plt.plot(ab,t2ccarray, 'k', )
#plt.plot (ab , t2cc , 'b--', label = 'Analytical')
plt.xlabel( '$Y_s$' )
plt.ylabel(' $T_2 ^{cc} $')
plt.title (' $T_2^ {cc}$')
plt.grid()
plt.legend()

plt.figure(7)

plt.plot(ab,t2csarray, 'b', )
#plt.plot (ab , t2cs , 'b--', label = 'Analytical')
plt.xlabel( '$Y_s$' )
plt.ylabel(' $T_2 ^{cs} $')
plt.title ('Quadrupole coefficient $T_2^ {cs}$')
plt.grid()
plt.legend()

plt.figure(8)

plt.plot(ab,t2scarray, 'b', )
#plt.plot (ab , t2sc , 'b--', label = 'Analytical')
plt.xlabel( '$Y_s$' )
plt.ylabel(' $T_2 ^{sc} $')
plt.title ('Quadrupole coefficient $T_2^ {sc}$')
plt.grid()
plt.legend()

#plt.subplot (3,3,9)

#plt.plot(ab,t2ssarray, 'r--', label = 'Matrix')
#plt.plot (ab , t2ss , 'b--', label = 'Analytical')
#plt.xlabel( '$Y_s$' )
#plt.ylabel(' $T_2 ^{ss} $')
#plt.title ('Quadrupole coefficient $T_2^ {ss}$')
#plt.grid()
#plt.legend()

plt.figure(10)

plt.plot(ab,crossC, 'b', )
#plt.xlim([0, 350])
#plt.plot (ab , t2ss , 'b--', label = 'Analytical')
plt.xlabel( '$Y_s$' )
plt.ylabel(' cross section  ')
plt.title ('Cross Section for  compressional  ')
plt.grid()
#plt.legend()
plt.figure(11)
plt.plot(ab.real,t0ccreal, 'b', )
#plt.plot(ab.imag,t0ccimag, 'b', )
plt.xlabel( 'Real $Y_s$' )
plt.ylabel(' Reall $T_0^{cc}$  ')
plt.title ('Real part of Monopole $T_0^{cc}$  ')
plt.grid()

plt.figure(12)
plt.plot(ab.imag,t0ccimag, 'b', )
plt.xlabel( 'Imaginary $Y_s$' )
plt.ylabel(' Imaginary $T_0^{cc}$  ')
plt.title ('Imaginary part of Monopole $T_0^{cc}$  ')
plt.grid()

plt.figure(13)
plt.plot(ab.real,t1ccreal, 'b', )
plt.xlabel( 'Real $Y_s$' )
plt.ylabel(' Real $T_1^{cc}$  ')
plt.title ('Real part of Dipole $T_1^{cc}$  ')
plt.grid()

plt.figure(14)
plt.plot(ab.imag,t1ccimag, 'b', )
plt.xlabel( 'Imaginary $Y_s$' )
plt.ylabel(' Imaginary $T_1^{cc}$  ')
plt.title ('Imaginary part of Dipole $T_1^{cc}$  ')
plt.grid()

plt.figure(15)
plt.plot(ab.real,t1csreal, 'b', )
plt.xlabel( 'Real $Y_s$' )
plt.ylabel(' Real $T_1^{cs}$  ')
plt.title ('Real part of Dipole $T_1^{cs}$  ')
plt.grid()

plt.figure(16)
plt.plot(ab.imag,t1csimag, 'b', )
plt.xlabel( 'Imaginary $Y_s$' )
plt.ylabel(' Imaginary $T_1^{cs}$  ')
plt.title ('Imaginary part of Dipole $T_1^{cs}$  ')
plt.grid()

plt.figure(17)
plt.plot(ab.real,t2ccreal, 'b', )
plt.xlabel( 'Real $Y_s$' )
plt.ylabel(' Real $T_2^{cc}$  ')
plt.title ('Real part of quadrupole $T_2^{cc}$  ')
plt.grid()

plt.figure(18)
plt.plot(ab.imag,t2ccimag, 'b', )
plt.xlabel( 'Imaginary $Y_s$' )
plt.ylabel(' Imaginary $T_2^{cc}$  ')
plt.title ('Imaginary part of Quadrupole $T_2^{cc}$  ')
plt.grid()

plt.figure(19)
plt.plot(ab.real,t2csreal, 'b', )
plt.xlabel( 'Real $Y_s$' )
plt.ylabel(' Real $T_2^{cs}$  ')
plt.title ('Real part of Quadrupole $T_2^{cs}$  ')
plt.grid()

plt.figure(20)
plt.plot(ab.imag,t2csimag, 'b', )
plt.xlabel( 'Imaginary $Y_s$' )
plt.ylabel(' Imaginary $T_2^{cs}$  ')
plt.title ('Imaginary part of Quadrupole $T_2^{cs}$  ')
plt.grid()




    
    
    



    
