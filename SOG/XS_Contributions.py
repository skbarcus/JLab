import numpy as np
import matplotlib.pyplot as plt
#from scipy.interpolate import make_interp_spline
#from scipy.interpolate import interp1d
# Import curve fitting package from scipy
from scipy.optimize import curve_fit
#Import polynomial fitting from numpy.
from numpy.polynomial import Polynomial
#import sympy as sym
from scipy.misc import derivative
import scipy.integrate as integrate
from scipy.integrate import quad
from scipy.stats import norm

pi = 3.141592654
deg2rad = pi/180.0
GeV2fm = 1./0.0389                  #Convert Q^2 units from GeV^2 to fm^-2.
hbar = 6.582*np.power(10.0,-16.0)   #hbar in [eV*s].
C = 299792458.0                     #Speed of light [m/s]. 
e = 1.60217662E-19                  #Electron charge C.
alpha = 0.0072973525664             #1.0/137.0 Fine structure constant.
muHe3 = -2.1275*(3.0/2.0)           #Diens has this 3/2 factor for some reason, but it fits the data much better.  #2*2.793-1.913 is too naive.
ngaus = 12                          #Number of Gaussians used to fit data.
Z = 2.                              #Atomic number He3.
A = 3.                              #Mass number He3.
MtHe3 = 3.0160293*0.9315            #Mass of He3 in GeV.
gamma = 0.8*np.power(2.0/3.0,0.5)   #Gaussian width [fm] from Amroun gamma*sqrt(3/2) = 0.8 fm.
theta = 0.#21.04;
theta_cor = 0.                 #Theta that corrects for the Q^2eff adjustment. Basically when we plot the XS and FFs the Q2[0] is really Q^2eff if we don't use this theta_cor. This variable is for the slightly smaller theta representing the real scattering angle.
E0 = 1                     #3.356 Initial e- energy GeV.#2.216
Ef = 0. 

#My 3He thesis values.
R = (0.3, 0.7, 0.9, 1.1, 1.5, 1.6, 2.2, 2.7, 3.3, 4.2, 4.3, 4.8)
Qich = (0.0996392,0.214304,0.0199385,0.195676,0.0785533,0.167223,0.126926,0.0549379,0.0401401,0.0100803,0.0007217,4.98962e-12)
Qim = (0.159649,0.0316168,0.277843,0.0364955,0.0329718,0.233469,0.117059,0.0581085,0.0485212,1.77602e-12,0.0240927,8.94934e-12)

#Define charge form factor function.
def Fch(Q2eff,Qich,R):
    sumFch_ff = 0
    Fch_ff = 0
    for i in range(ngaus):
        sumFch_ff = (Qich[i]/(1.0+2.0*np.power(R[i],2.0)/np.power(gamma,2.0))) * ( np.cos(np.power(Q2eff,0.5)*R[i]) + (2.0*np.power(R[i],2.0)/np.power(gamma,2.0)) * (np.sin(np.power(Q2eff,0.5)*R[i])/(np.power(Q2eff,0.5)*R[i])) )
        Fch_ff = Fch_ff + sumFch_ff
    Fch_ff =  Fch_ff * np.exp(-0.25*Q2eff*np.power(gamma,2.0))
    return Fch_ff

#Define magnetic form factor function.
def Fm(Q2eff,Qim,R):
    sumFm_ff = 0
    Fm_ff = 0
    for i in range(ngaus):
        sumFm_ff = (Qim[i]/(1.0+2.0*np.power(R[i],2.0)/np.power(gamma,2.0))) * ( np.cos(np.power(Q2eff,0.5)*R[i]) + (2.0*np.power(R[i],2.0)/np.power(gamma,2.0)) * (np.sin(np.power(Q2eff,0.5)*R[i])/(np.power(Q2eff,0.5)*R[i])) )
        Fm_ff = Fm_ff + sumFm_ff
    Fm_ff =  Fm_ff * np.exp(-0.25*Q2eff*np.power(gamma,2.0))
    return Fm_ff

def XS(Q2eff,E0):

    theta = 2*np.arcsin(  np.power( (1/(4*np.power(E0,2.)*GeV2fm/Q2eff-2*E0/MtHe3)) , 0.5 )  )
    theta_cor = 2*np.arcsin(  np.power( (1/((4*np.power(E0,2.)*GeV2fm/   np.power(  np.power(Q2eff,0.5)/(1+(1.5*2*alpha)/(E0*np.power(GeV2fm,0.5)*1.12*np.power(3.,1./3.)))  ,2.)   )-(2.*E0/MtHe3))) , 0.5 )  )
    Ef = E0/(1.0+2.0*E0*np.power(np.sin(theta/2.0),2.0)/MtHe3)
    W = E0 - Ef
    q2_3 = abs(  np.power(W,2.0)*GeV2fm - Q2eff  )        #Convert w^2 from GeV^2 to fm^-2 to match Q2. [fm^-2]
    eta = 1.0 + Q2eff/(4.0*np.power(MtHe3,2.0)*GeV2fm)       #Make sure Mt^2 is converted from GeV^2 to fm^-2 to match Q^2. 

    #Calculate Mott XS.
    mottxs = (  (np.power(Z,2.)*(Ef/E0)) * (np.power(alpha,2.0)/(4.0*np.power(E0,2.0)*np.power(np.sin(theta/2.0),4.0)))*np.power(np.cos(theta/2.0),2.0)  ) * 1.0/25.7    #Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
    
    #Calculate XS from FFs.
    xs = mottxs * (1./eta) * ( (Q2eff/q2_3)*np.power(Fch(Q2eff,Qich,R),2.) + (np.power(muHe3,2.0)*Q2eff/(2*np.power(MtHe3,2)*GeV2fm))*(0.5*Q2eff/q2_3 + np.power(np.tan(theta/2),2))*np.power(Fm(Q2eff,Qim,R),2.) )
    xs_ch = mottxs * (1./eta) * ( (Q2eff/q2_3)*np.power(Fch(Q2eff,Qich,R),2.) )
    xs_m = mottxs * (1./eta) * ( (np.power(muHe3,2.0)*Q2eff/(2*np.power(MtHe3,2)*GeV2fm))*(0.5*Q2eff/q2_3 + np.power(np.tan(theta/2),2))*np.power(Fm(Q2eff,Qim,R),2.) )

    ratio_ch = xs_ch / xs
    ratio_m = xs_m / xs

    return xs, xs_ch, xs_m, ratio_ch, ratio_m; 

print(XS(10,E0))

#Plot the XS along with the contributions from the charge and magnetic parts.
fig, ax = plt.subplots(figsize=(12,6))
ax.set_title('$^3$He Cross Section at {:.3f} GeV'.format(E0),fontsize=20)
ax.set_ylabel(r'$\frac{d\sigma}{d\Omega}$ (fm$^2$/sr)',fontsize=16)
ax.set_xlabel('$Q^2$ (fm$^{-2}$)',fontsize=16)
ax.set_yscale('log')

min = 0
max = 60
plt.xticks(np.arange(min, max+1, 2.0))

#Define Q2eff range to plot.
Q2eff = np.linspace(0.00001,60,600)

plt.plot(Q2eff, XS(Q2eff,E0)[0], color='red', alpha=1.)
plt.plot(Q2eff, XS(Q2eff,E0)[1], color='blue', alpha=1.)
plt.plot(Q2eff, XS(Q2eff,E0)[2], color='green', alpha=1.)

plt.show()

#Plot the relative contributions of the charge and magnetic FFs.
fig, ax = plt.subplots(figsize=(12,6))
ax.set_title('$^3$He Cross Section Charge and Magnetic Contributions at {:.3f} GeV'.format(E0),fontsize=20)
ax.set_ylabel('Relative Contribution to Cross Section',fontsize=16)
ax.set_xlabel('$Q^2$ (fm$^{-2}$)',fontsize=16)
#ax.set_yscale('log')

min = 0
max = 60
plt.xticks(np.arange(min, max+1, 2.0))

#Define Q2eff range to plot.
Q2eff = np.linspace(0.00001,60,600)

plt.plot(Q2eff, XS(Q2eff,E0)[3], color='blue', alpha=1.)
plt.plot(Q2eff, XS(Q2eff,E0)[4], color='green', alpha=1.)

plt.show()
