import numpy as np
import scipy.special as sc

import matplotlib.pyplot as plt

T0ccComb = []
ysComb = []
T0ccReal = []
ysReal = []
ysImag = []
T0ccImag = []
a = 250*(10**-6) #radius of the particle
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

for f in range(10000, 10006, 1):
    
    kc = 2*pi*(f/vo) + 1j*att_o*f**2
    kcp = 2*pi*f/vi + 1j*att_i*f**2
    ks = (1+1j)*(np.sqrt(rho_o*2*pi*f/(2*mu_pp)))
    ksp =(np.sqrt(rho_i))*(2*pi*f)/(np.sqrt(mu0_s))


    yc = kc*a
    ycp = kcp*a
    ys = ks*a
    ysp = ksp*a
   

    j0yc =  sc.spherical_jn(0, yc,) #j0(yc)
   
    j0ycp =  sc.spherical_jn(0, ycp,) #j0(yc')
    
    j0pyc = sc.spherical_jn(0, yc, [1]) #j'0(yc)
   

    j0pycp = sc.spherical_jn(0, ycp, [1]) #j'0(yc')
    

    h0yc = sc.spherical_jn(0,yc) +(sc.spherical_yn(0,yc))*1j  #h0(yc)
    
    h0pyc = sc.spherical_jn(0, yc, [1])+ 1j*sc.spherical_yn(0, yc, [1]) #h'0(yc)
    




    # Calculating matrix component:

  
    # M11  = yc* h0'(yc) 
    M11 = yc*h0pyc
    # print ("M11  :" , M11)

    # M12  = - yc' j0'(yc') 
    M12 = -(ycp*j0pycp)
    # print ("M12  :" , M12)

    #M21   = 1/ys^2 [(ys^2 h0(yc)) + 4 yc h0'(yc)]
    M21 = ( (ys**2 * h0yc) + (4*yc * h0pyc)) * (1/ys**2)
    # print ("M21  :" , M21)

    # M22   = -rho_cap/ys'^2 [(ys'^2 j0(yc') + 4 yc' j0'(yc')]
    M22 = -((ysp**2 *j0ycp) + (4*ycp*j0pycp))*(rho_cap /ysp**2)
    # print ("M22  :" , M22)

    # C1  = - yc j0'(yc)
    C1 = -(yc*j0pyc)
    # print ("C1  :" , C1)

    # C2  = -1/ys^2 [(ys^2 j0(yc) + 4 yc j0'(yc)]
    C2 = -((ys**2 *j0yc)+ (4*yc *j0pyc))*(1/ys**2)
    # print ("C2  :" , C2)

   

    M = np.array([
        [M11,M12],
        [M21,M22]])
#     print (np.shape(M))
   
    
   

    C = np.array([[C1],[C2]])
   



    inverse=np.linalg.inv(M)
    X = np.dot(inverse,C)
   

    
    T0cc = X[0]

    T0ccReal.append(T0cc.real)
    ysReal.append(ys.real)
    T0ccImag.append (T0cc.imag)
    ysImag.append (ys.imag)
    
    T0ccComb.append(T0cc)
    ysComb.append (ys)
    
print (ysComb)
print (T0ccComb)
    
    
   


# plt.figure(figsize=(12,8))

# plt.subplot (2,2,1)
# plt.plot(ysComb, T0ccComb, 'r--')
# plt.xlabel( ' $Y_s$' )
# plt.ylabel('$T_0 ^{cc} $')
# plt.title (' $T_0^ {cc}$')
# plt.grid()

# plt.subplot (2,2,2)
# plt.plot(ysImag, T0ccImag, 'r--')
# plt.xlabel('Imag $Y_s$')
# plt.ylabel('Imag $T_0 ^{cc} $')
# plt.title ('Imaginary Part of $T_0^ {cc}$')
# plt.grid()



# plt.tight_layout()
# plt.show()


   