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
rad2deg = 180.0/pi
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
max_He3_Q2_val = 4*np.power(E0,2) / (1 + 2*E0/MtHe3) * GeV2fm #Calculate the maximum attainable Q2 value for 3He for the given beam energy.

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

def Q2_2_theta(Q2eff):
    #print('Q2 =',Q2eff,' Theta =', 2*np.arcsin(  np.power( (1/(4*np.power(E0,2.)*GeV2fm/Q2eff-2*E0/MtHe3)) , 0.5 )  )  * rad2deg)
    theta =  2*np.arcsin(  np.power( (1/(4*np.power(E0,2.)*GeV2fm/Q2eff-2*E0/MtHe3)) , 0.5 )  )  * rad2deg
    return ["%.0f" % z for z in theta]

def Q2_2_theta_vals(Q2eff):
    #print('Q2 =',Q2eff,' Theta =', 2*np.arcsin(  np.power( (1/(4*np.power(E0,2.)*GeV2fm/Q2eff-2*E0/MtHe3)) , 0.5 )  )  * rad2deg)
    theta =  2*np.arcsin(  np.power( (1/(4*np.power(E0,2.)*GeV2fm/Q2eff-2*E0/MtHe3)) , 0.5 )  )  * rad2deg
    return theta

def theta_2_Q2(theta):
    #print('theta =',theta,' Q2 =',4*E0*E0*np.power(np.sin(theta*deg2rad/2.0),2.0) * GeV2fm / (1.0+2.0*E0*np.power(np.sin(theta*deg2rad/2.0),2.0)/MtHe3))
    return 4*E0*E0*np.power(np.sin(theta*deg2rad/2.0),2.0) * GeV2fm / (1.0+2.0*E0*np.power(np.sin(theta*deg2rad/2.0),2.0)/MtHe3)

###########################################################
#Plot the XS and the charge and magnetic FF contributions.#
###########################################################

#Define Q2eff range to plot.
minq2 = 0.00001 #Avoid poles at zero.
if max_He3_Q2_val < 60:
    maxq2 = max_He3_Q2_val-0.1
else:
    maxq2 = 60
nsteps = 1000
dstep = 5
Q2eff = np.linspace(minq2,maxq2,nsteps)

print("Q2_2_theta_vals(0.959)",Q2_2_theta_vals(0.959))

#Plot the XS along with the contributions from the charge and magnetic parts.
fig, ax1 = plt.subplots(figsize=(12,7))
ax1.set_title('$^3$He Cross Section', fontsize=20)
#ax1.set_ylabel(r'$\frac{d\sigma}{d\Omega}$ (fm$^2$/sr)',fontsize=16)
ax1.set_ylabel(r'$d\sigma/d\Omega$ (fm$^2$/sr)',fontsize=16)
ax1.set_xlabel('$Q^2$ (fm$^{-2}$)',fontsize=16)
ax1.set_yscale('log')
#Define secondary X-axis twinned off of first one's Y-axis.
ax2 = ax1.twiny()
ax2.set_xlabel(r"Scattering Angle ($^\circ$)",fontsize=16)
ax2.xaxis.set_label_coords(.48, .95)

ax1.plot(Q2eff, XS(Q2eff,E0)[0], color='red', label='$^3$He Cross Section',alpha=1.,linewidth=2,zorder=2.5)
ax1.plot(Q2eff, XS(Q2eff,E0)[1], color='blue', label='Electric Form Factor Cross Section Contribution', alpha=1.,linewidth=2.5,zorder=1)
ax1.plot(Q2eff, XS(Q2eff,E0)[2], color='green', label='Magnetic Form Factor Cross Section Contribution', alpha=1.,linewidth=2.5,zorder=0)

new_tick_locations = np.arange(minq2,maxq2+5,step=dstep)
#new_tick_locations = np.arange(minq2,maxq2+5,step=dstep) #1GeV nice scattering angle ticks.
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(Q2_2_theta(new_tick_locations))

ax2.plot([], [], ' ', label='Initial Electron Energy = {:.3f} GeV'.format(E0))
ax2.plot([], [], ' ', label='Maximum $Q^2$ = {:.3f} fm$^2$'.format(max_He3_Q2_val))
ax1.legend(loc='center right',fontsize=12)
ax2.legend(loc='upper right',fontsize=12)
plt.show()

#################################################################
#Plot the relative contributions of the charge and magnetic FFs.#
#################################################################
fig, ax1 = plt.subplots(figsize=(12,7))
ax1.set_title('$^3$He Cross Section Charge and Magnetic Contributions',fontsize=20)
ax1.set_ylabel('Relative Contribution to Cross Section',fontsize=16)
ax1.set_xlabel('$Q^2$ (fm$^{-2}$)',fontsize=16)
#ax1.set_yscale('log')

#Define secondary X-axis twinned off of first one's Y-axis.
ax2 = ax1.twiny()
ax2.set_xlabel(r"Scattering Angle ($^\circ$)",fontsize=16)
ax2.xaxis.set_label_coords(.48, .95)

new_tick_locations = np.arange(minq2,maxq2,step=dstep)
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(Q2_2_theta(new_tick_locations))

plt.xticks(np.arange(minq2, maxq2, dstep))

ax1.plot(Q2eff, XS(Q2eff,E0)[3], color='blue', label='Electric Form Factor', alpha=1.)
ax1.plot(Q2eff, XS(Q2eff,E0)[4], color='green', label='Magnetic Form Factor', alpha=1.)

# Create empty plot with blank marker containing extra label.
ax2.plot([], [], ' ', label='Initial Electron Energy = {:.3f} GeV'.format(E0))
ax2.plot([], [], ' ', label='Maximum $Q^2$ = {:.3f} fm$^2$'.format(max_He3_Q2_val))
ax1.legend(loc='center left',title='Fraction of Cross Section Due to:',fontsize=12)
ax2.legend(loc='upper right',fontsize=12)
plt.show()
