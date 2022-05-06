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

#use_3He = 0                         #If plotting 3He data. 
#use_3H  = 1                         #
pi = 3.141592654
deg2rad = pi/180.0
rad2deg = 180.0/pi
GeV2fm = 1./0.0389                  #Convert Q^2 units from GeV^2 to fm^-2.
hbar = 6.582*np.power(10.0,-16.0)   #hbar in [eV*s].
C = 299792458.0                     #Speed of light [m/s]. 
e = 1.60217662E-19                  #Electron charge C.
alpha = 0.0072973525664             #1.0/137.0 Fine structure constant.
muHe3 = -2.1275*(3.0/2.0)           #Diens has this 3/2 factor for some reason, but it fits the data much better.  #2*2.793-1.913 is too naive.
muH3 = 2.9788*(3.0/1.0)             #Magnetic moment of trinucleon (H3 or He3). NIST: http://physics.nist.gov/cgi-bin/cuu/Results?search_for=magnet+moment   //MCEEP Code for H3 and He3 eleastic FFs has magnetic moments multiplied by 3.0/Z. I don't know why but it works. Maybe it's a factor of A/Z?
ngaus = 12                          #Number of Gaussians used to fit data.
Z = 2.                              #Atomic number He3.
A = 3.                              #Mass number He3.
MtHe3 = 3.0160293*0.9315            #Mass of He3 in GeV.
MtH3 = 3.0160492*0.9315;            #Mass of H3 in GeV.
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

#My 3H thesis values.
R_H3_thesis = (0.3,0.8,1.4,1.9,2.5,3.3,4.1,4.8)
Qich_H3_thesis = (0.151488,0.348372,0.29635,0.0978631,0.121983,0.0242654,0.049329,4.40751e-11)
Qim_H3_thesis = (0.190646,0.301416,0.318972,0.159433,0.173933,0.106361,0.0665564,0.0148866)

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
    
    #Calculate XS from FFs. Divide out Mott for pretty plots.
    xs = (1./eta) * ( (Q2eff/q2_3)*np.power(Fch(Q2eff,Qich,R),2.) + (np.power(muHe3,2.0)*Q2eff/(2*np.power(MtHe3,2)*GeV2fm))*(0.5*Q2eff/q2_3 + np.power(np.tan(theta/2),2))*np.power(Fm(Q2eff,Qim,R),2.) )
    xs_ch = (1./eta) * ( (Q2eff/q2_3)*np.power(Fch(Q2eff,Qich,R),2.) )
    xs_m = (1./eta) * ( (np.power(muHe3,2.0)*Q2eff/(2*np.power(MtHe3,2)*GeV2fm))*(0.5*Q2eff/q2_3 + np.power(np.tan(theta/2),2))*np.power(Fm(Q2eff,Qim,R),2.) )

    ratio_ch = xs_ch / xs
    ratio_m = xs_m / xs

    return xs, xs_ch, xs_m, ratio_ch, ratio_m, Q2eff, Ef, theta*180/pi; 

def XS_Data(E0,theta):

    #theta = 2*np.arcsin(  np.power( (1/(4*np.power(E0,2.)*GeV2fm/Q2eff-2*E0/MtHe3)) , 0.5 )  )
    #theta_cor = 2*np.arcsin(  np.power( (1/((4*np.power(E0,2.)*GeV2fm/   np.power(  np.power(Q2eff,0.5)/(1+(1.5*2*alpha)/(E0*np.power(GeV2fm,0.5)*1.12*np.power(3.,1./3.)))  ,2.)   )-(2.*E0/MtHe3))) , 0.5 )  )
    Ef = E0/(1.0+2.0*E0*np.power(np.sin(theta*deg2rad/2.0),2.0)/MtHe3)
    Q2 = 4.0*E0*Ef*np.power(np.sin(theta*deg2rad/2.0),2.0) * GeV2fm
    Q2eff = np.power( np.power(Q2,0.5) * (1.0+(1.5*Z*alpha)/(E0*np.power(GeV2fm,0.5)*1.12*np.power(A,1.0/3.0))) ,2.0)

    W = E0 - Ef
    q2_3 = abs(  np.power(W,2.0)*GeV2fm - Q2eff  )        #Convert w^2 from GeV^2 to fm^-2 to match Q2. [fm^-2]
    eta = 1.0 + Q2eff/(4.0*np.power(MtHe3,2.0)*GeV2fm)       #Make sure Mt^2 is converted from GeV^2 to fm^-2 to match Q^2. 

    tau = np.power(Q2eff,2)/(4*np.power(MtHe3,2))
    epsilon = np.power( 1+2*(1+tau)*np.power(np.tan(theta/2*deg2rad),2) ,-1)

    #Calculate Mott XS.
    mottxs = (  (np.power(Z,2.)*(Ef/E0)) * (np.power(alpha,2.0)/(4.0*np.power(E0,2.0)*np.power(np.sin(theta/2.0*deg2rad),4.0)))*np.power(np.cos(theta/2.0*deg2rad),2.0)  ) * 1.0/25.7    #Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
    
    #Calculate XS from FFs. Divide out Mott for pretty plots.
    xs = mottxs * (1./eta) * ( (Q2eff/q2_3)*np.power(Fch(Q2eff,Qich,R),2.) + (np.power(muHe3,2.0)*Q2eff/(2*np.power(MtHe3,2)*GeV2fm))*(0.5*Q2eff/q2_3 + np.power(np.tan(theta/2),2))*np.power(Fm(Q2eff,Qim,R),2.) )
    xs_ch = mottxs * (1./eta) * ( (Q2eff/q2_3)*np.power(Fch(Q2eff,Qich,R),2.) )
    xs_m = mottxs * (1./eta) * ( (np.power(muHe3,2.0)*Q2eff/(2*np.power(MtHe3,2)*GeV2fm))*(0.5*Q2eff/q2_3 + np.power(np.tan(theta/2),2))*np.power(Fm(Q2eff,Qim,R),2.) )

    ratio_ch = xs_ch / xs
    ratio_m = xs_m / xs

    return xs, xs_ch, xs_m, ratio_ch, ratio_m, Q2eff, Ef, theta, epsilon; 


def XS_Data_H3(E0,theta):
    ngaus = 8      #Number of Gaussians used to fit 3H data.
    Z = 1          #Atomic number H3.

    #My 3H thesis values.
    R = (0.3,0.8,1.4,1.9,2.5,3.3,4.1,4.8)
    Qich = (0.151488,0.348372,0.29635,0.0978631,0.121983,0.0242654,0.049329,4.40751e-11)
    Qim = (0.190646,0.301416,0.318972,0.159433,0.173933,0.106361,0.0665564,0.0148866)

    Ef = E0/(1.0+2.0*E0*np.power(np.sin(theta*deg2rad/2.0),2.0)/MtH3)
    Q2 = 4.0*E0*Ef*np.power(np.sin(theta*deg2rad/2.0),2.0) * GeV2fm
    Q2eff = np.power( np.power(Q2,0.5) * (1.0+(1.5*Z*alpha)/(E0*np.power(GeV2fm,0.5)*1.12*np.power(A,1.0/3.0))) ,2.0)

    W = E0 - Ef
    q2_3 = abs(  np.power(W,2.0)*GeV2fm - Q2eff  )        #Convert w^2 from GeV^2 to fm^-2 to match Q2. [fm^-2]
    eta = 1.0 + Q2eff/(4.0*np.power(MtH3,2.0)*GeV2fm)       #Make sure Mt^2 is converted from GeV^2 to fm^-2 to match Q^2. 

    tau = np.power(Q2eff,2)/(4*np.power(MtH3,2))
    epsilon = np.power( 1+2*(1+tau)*np.power(np.tan(theta/2*deg2rad),2) ,-1)

    #Calculate Mott XS.
    mottxs = (  (np.power(Z,2.)*(Ef/E0)) * (np.power(alpha,2.0)/(4.0*np.power(E0,2.0)*np.power(np.sin(theta/2.0*deg2rad),4.0)))*np.power(np.cos(theta/2.0*deg2rad),2.0)  ) * 1.0/25.7    #Convert GeV^-2 to fm^2 by multiplying by 1/25.7.
    
    #Calculate XS from FFs. Divide out Mott for pretty plots.
    xs = mottxs * (1./eta) * ( (Q2eff/q2_3)*np.power(Fch(Q2eff,Qich,R),2.) + (np.power(muH3,2.0)*Q2eff/(2*np.power(MtH3,2)*GeV2fm))*(0.5*Q2eff/q2_3 + np.power(np.tan(theta/2),2))*np.power(Fm(Q2eff,Qim,R),2.) )
    xs_ch = mottxs * (1./eta) * ( (Q2eff/q2_3)*np.power(Fch(Q2eff,Qich,R),2.) )
    xs_m = mottxs * (1./eta) * ( (np.power(muH3,2.0)*Q2eff/(2*np.power(MtH3,2)*GeV2fm))*(0.5*Q2eff/q2_3 + np.power(np.tan(theta/2),2))*np.power(Fm(Q2eff,Qim,R),2.) )

    ratio_ch = xs_ch / xs
    ratio_m = xs_m / xs

    return xs, xs_ch, xs_m, ratio_ch, ratio_m, Q2eff, Ef, theta, epsilon; 

#print('XS predicted for my 3He data point (Q2,E0) =',XS(34.19,3.356))

#print('XS predicted for my 3He data point (E0,theta) =',XS_Data(3.356,20.51))

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

#print("Q2_2_theta_vals(0.959)",Q2_2_theta_vals(0.959))

#Plot the XS along with the contributions from the charge and magnetic parts.
fig, ax1 = plt.subplots(figsize=(12,7))
ax1.set_title('$^3$He Cross Section', fontsize=20)
#ax1.set_ylabel(r'$\frac{d\sigma}{d\Omega}$ (fm$^2$/sr)',fontsize=16)
ax1.set_ylabel(r'$\left( d\sigma/d\Omega \right) \; / \; \left( d\sigma/d\Omega \right)_{Mott}$',fontsize=16)
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
ax1.legend(loc=(0.5,0.65),fontsize=12)
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

###############################################################################################
#Plot the charge and magnetic fractional contributions to the total XS for the 3He world data.#
###############################################################################################

#Read in data.
with open('/home/skbarcus/JLab/SOG/3He_Data.txt') as f:
#with open('/home/skbarcus/JLab/SOG/3H_Data_Thesis_Datasets.txt') as f:
#with open('/home/skbarcus/JLab/SOG/3H_Data_Thesis.txt') as f:
    lines = f.readlines()

#Remove lines with column labels.
del lines[0]
del lines[0]

#Create arrays.
Raw_Data_He3 = []          #Array to store raw data.

#Read each line and split by the spaces and store as array.
#['Energy (GeV)', 'Theta (Degrees)', 'Sigma Experimental', 'Uncertainties', 'Dataset']
for line in lines:
    event = line.split()
    Raw_Data_He3.append(event)

#Turn data into numpy array and swap char to float.
Raw_Data_He3 = np.array(Raw_Data_He3)
Raw_Data_He3 = np.array(Raw_Data_He3.astype('float'))

#Examine data shape and check output.
#print('Raw_Data_He3.shape = ',Raw_Data_He3.shape)
#print('Raw_Data_He3[0] = ',Raw_Data_He3[0])

XS_Data_Q2eff_He3 = [[],[],[],[],[],[],[],[]]
XS_Data_ch_frac_He3 = [[],[],[],[],[],[],[],[]]
XS_Data_mag_frac_He3 = [[],[],[],[],[],[],[],[]]
XS_Data_Residual_He3 = [[],[],[],[],[],[],[],[]]
XS_Data_Epsilon_He3 = [[],[],[],[],[],[],[],[]]

for i in range(0,len(Raw_Data_He3)):
    #Calculate the residual per data point using the SOG fit result.
    res = (Raw_Data_He3[i][2] - XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[0] ) /  XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[0]
    #print('XSexp =',Raw_Data_He3[i][2],'   XSfit =',XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[0],'   Residual =',res)
    if Raw_Data_He3[i][4]==1:#Amroun 1994
        XS_Data_Q2eff_He3[0].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[5])
        XS_Data_ch_frac_He3[0].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[3])
        XS_Data_mag_frac_He3[0].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[4])
        XS_Data_Residual_He3[0].append(res)
        XS_Data_Epsilon_He3[0].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[8])
    if Raw_Data_He3[i][4]==2:#Collard 1965
        XS_Data_Q2eff_He3[1].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[5])
        XS_Data_ch_frac_He3[1].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[3])
        XS_Data_mag_frac_He3[1].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[4])
        XS_Data_Residual_He3[1].append(res)
        XS_Data_Epsilon_He3[1].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[8])
        #if XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[4]>0.75:
            #print(Raw_Data_He3[i])
    if Raw_Data_He3[i][4]==3:#Szalata 1977
        XS_Data_Q2eff_He3[2].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[5])
        XS_Data_ch_frac_He3[2].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[3])
        XS_Data_mag_frac_He3[2].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[4])
        XS_Data_Residual_He3[2].append(res)
        XS_Data_Epsilon_He3[2].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[8])
    if Raw_Data_He3[i][4]==4:#Dunn 1983
        XS_Data_Q2eff_He3[3].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[5])
        XS_Data_ch_frac_He3[3].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[3])
        XS_Data_mag_frac_He3[3].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[4])
        XS_Data_Residual_He3[3].append(res)
        XS_Data_Epsilon_He3[3].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[8])
    if Raw_Data_He3[i][4]==5:#Camsonne 2016
        XS_Data_Q2eff_He3[4].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[5])
        XS_Data_ch_frac_He3[4].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[3])
        XS_Data_mag_frac_He3[4].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[4])
        XS_Data_Residual_He3[4].append(res)
        XS_Data_Epsilon_He3[4].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[8])
    if Raw_Data_He3[i][4]==6:#Nakagawa 2001
        XS_Data_Q2eff_He3[5].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[5])
        XS_Data_ch_frac_He3[5].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[3])
        XS_Data_mag_frac_He3[5].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[4])
        XS_Data_Residual_He3[5].append(res)
        XS_Data_Epsilon_He3[5].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[8])
    if Raw_Data_He3[i][4]==7:#Barcus 2019
        XS_Data_Q2eff_He3[6].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[5])
        XS_Data_ch_frac_He3[6].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[3])
        XS_Data_mag_frac_He3[6].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[4])
        XS_Data_Residual_He3[6].append(res)
        XS_Data_Epsilon_He3[6].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[8])
    if Raw_Data_He3[i][4]==8:#Arnold 1978
        XS_Data_Q2eff_He3[7].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[5])
        XS_Data_ch_frac_He3[7].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[3])
        XS_Data_mag_frac_He3[7].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[4])
        XS_Data_Residual_He3[7].append(res)
        XS_Data_Epsilon_He3[7].append(XS_Data(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[8])

XS_Data_Q2eff_He3 = np.array(XS_Data_Q2eff_He3)
XS_Data_ch_frac_He3= np.array(XS_Data_ch_frac_He3)
XS_Data_mag_frac_He3= np.array(XS_Data_mag_frac_He3)

#print(XS_Data(Raw_Data_He3[0],Raw_Data_He3[1])[5])
#print(XS_Data_Q2eff_He3)
print('Raw_Data_He3.shape =',Raw_Data_He3.shape)
print('XS_Data_Q2eff_He3.shape =',XS_Data_Q2eff_He3.shape)

#Plot the world data by charge and magnetic XS fractions.
#Plot the charge contributions.
fig, ax = plt.subplots(figsize=(12,6))
ax.set_ylabel('Fractional Charge Contribution to Cross Section',fontsize=18)
ax.set_xlabel(r'$Q^2$ (fm$^{-2}$)',fontsize=18)
ax.set_title(r'Charge Contribution to the $^3$He Cross Section',fontsize=20)

#print(XS_Data_Q2eff_He3[0],XS_Data_ch_frac_He3[0])
plt.plot(XS_Data_Q2eff_He3[1],XS_Data_ch_frac_He3[1],'ro',label='Collard 1965')#Collard 1965
plt.plot(XS_Data_Q2eff_He3[2],XS_Data_ch_frac_He3[2],'go',color='saddlebrown',label='Szalata 1977')#Szalata 1977
plt.plot(XS_Data_Q2eff_He3[7],XS_Data_ch_frac_He3[7],'yo',label='Arnold 1978')#Arnold 1978
plt.plot(XS_Data_Q2eff_He3[3],XS_Data_ch_frac_He3[3],'co',label='Dunn 1983')#Dunn 1983
plt.plot(XS_Data_Q2eff_He3[0],XS_Data_ch_frac_He3[0],'bo',label='Amroun 1994')#Amroun 1994
plt.plot(XS_Data_Q2eff_He3[5],XS_Data_ch_frac_He3[5],'o',color='g',label='Nakagawa 2001')#Nakagawa 2001
plt.plot(XS_Data_Q2eff_He3[4],XS_Data_ch_frac_He3[4],'mo',label='Camsonne 2016')#Camsonne 2016
plt.plot(XS_Data_Q2eff_He3[6],XS_Data_ch_frac_He3[6],'o',color='darkorange',label='Barcus 2019')#Barcus 2019

ax.legend(loc='best',fontsize=16,fancybox=True,framealpha=0.3)
plt.show()

#Plot the magnetic contributions.
fig, ax = plt.subplots(figsize=(12,6))
ax.set_ylabel('Fractional Magnetic Contribution to Cross Section',fontsize=18)
ax.set_xlabel(r'$Q^2$ (fm$^{-2}$)',fontsize=18)
ax.set_title(r'Magnetic Contribution to the $^3$He Cross Section',fontsize=20)

#print(XS_Data_Q2eff_He3[0],XS_Data_ch_frac_He3[0])
plt.plot(XS_Data_Q2eff_He3[1],XS_Data_mag_frac_He3[1],'ro',label='Collard 1965')#Collard 1965
plt.plot(XS_Data_Q2eff_He3[2],XS_Data_mag_frac_He3[2],'o',color='saddlebrown',label='Szalata 1977')#Szalata 1977
plt.plot(XS_Data_Q2eff_He3[7],XS_Data_mag_frac_He3[7],'yo',label='Arnold 1978')#Arnold 1978
plt.plot(XS_Data_Q2eff_He3[3],XS_Data_mag_frac_He3[3],'co',label='Dunn 1983')#Dunn 1983
plt.plot(XS_Data_Q2eff_He3[0],XS_Data_mag_frac_He3[0],'bo',label='Amroun 1994')#Amroun 1994
plt.plot(XS_Data_Q2eff_He3[5],XS_Data_mag_frac_He3[5],'o',color='g',label='Nakagawa 2001')#Nakagawa 2001
plt.plot(XS_Data_Q2eff_He3[4],XS_Data_mag_frac_He3[4],'mo',label='Camsonne 2016')#Camsonne 2016
plt.plot(XS_Data_Q2eff_He3[6],XS_Data_mag_frac_He3[6],'o',color='darkorange',label='Barcus 2019')#Barcus 2019

ax.legend(loc='best',fontsize=16,fancybox=True,framealpha=0.3)
plt.show()


#Plot Q2 vs. epsilon for each data point.
fig, ax = plt.subplots(figsize=(12,6))
ax.set_ylabel('Epsilon',fontsize=18)
ax.set_xlabel(r'$Q^2$ (fm$^{-2}$)',fontsize=18)
ax.set_title(r'World Data $Q^2$ vs. Epsilon',fontsize=20)

outer_rad = 80 #Marker size = (diameter in points)^2 
"""
size0 = [outer_rad*n for n in XS_Data_ch_frac_He3[0]]
size1 = [outer_rad*n for n in XS_Data_ch_frac_He3[1]]
size2 = [outer_rad*n for n in XS_Data_ch_frac_He3[2]]
size3 = [outer_rad*n for n in XS_Data_ch_frac_He3[3]]
size4 = [outer_rad*n for n in XS_Data_ch_frac_He3[4]]
size5 = [outer_rad*n for n in XS_Data_ch_frac_He3[5]]
size6 = [outer_rad*n for n in XS_Data_ch_frac_He3[6]]
size7 = [outer_rad*n for n in XS_Data_ch_frac_He3[7]]
"""
size0 = [outer_rad*n*pi/4 for n in XS_Data_mag_frac_He3[0]]
size1 = [outer_rad*n*pi/4 for n in XS_Data_mag_frac_He3[1]]
size2 = [outer_rad*n*pi/4 for n in XS_Data_mag_frac_He3[2]]
size3 = [outer_rad*n*pi/4 for n in XS_Data_mag_frac_He3[3]]
size4 = [outer_rad*n*pi/4 for n in XS_Data_mag_frac_He3[4]]
size5 = [outer_rad*n*pi/4 for n in XS_Data_mag_frac_He3[5]]
size6 = [outer_rad*n*pi/4 for n in XS_Data_mag_frac_He3[6]]
size7 = [outer_rad*n*pi/4 for n in XS_Data_mag_frac_He3[7]]

#print(XS_Data_Q2eff_He3[0],XS_Data_ch_frac_He3[0])
plt.scatter(XS_Data_Q2eff_He3[1],XS_Data_Epsilon_He3[1],edgecolors='r',facecolors='none',s=outer_rad,label='Collard 1965')#Collard 1965
plt.scatter(XS_Data_Q2eff_He3[2],XS_Data_Epsilon_He3[2],edgecolors='saddlebrown',facecolors='none',s=outer_rad,label='Szalata 1977')#Szalata 1977
plt.scatter(XS_Data_Q2eff_He3[7],XS_Data_Epsilon_He3[7],edgecolors='y',facecolors='none',s=outer_rad,label='Arnold 1978')#Arnold 1978
plt.scatter(XS_Data_Q2eff_He3[3],XS_Data_Epsilon_He3[3],edgecolors='c',facecolors='none',s=outer_rad,label='Dunn 1983')#Dunn 1983
plt.scatter(XS_Data_Q2eff_He3[0],XS_Data_Epsilon_He3[0],edgecolors='b',facecolors='none',s=outer_rad,label='Amroun 1994')#Amroun 1994
plt.scatter(XS_Data_Q2eff_He3[5],XS_Data_Epsilon_He3[5],edgecolors='g',facecolors='none',s=outer_rad,label='Nakagawa 2001')#Nakagawa 2001
plt.scatter(XS_Data_Q2eff_He3[4],XS_Data_Epsilon_He3[4],edgecolors='m',facecolors='none',s=outer_rad,label='Camsonne 2016')#Camsonne 2016
plt.scatter(XS_Data_Q2eff_He3[6],XS_Data_Epsilon_He3[6],edgecolors='darkorange',facecolors='none',s=outer_rad,label='Barcus 2019')#Barcus 2019

plt.scatter(XS_Data_Q2eff_He3[1],XS_Data_Epsilon_He3[1],color='r',s=size1)#Collard 1965
plt.scatter(XS_Data_Q2eff_He3[2],XS_Data_Epsilon_He3[2],color='saddlebrown',s=size2)#Szalata 1977
plt.scatter(XS_Data_Q2eff_He3[7],XS_Data_Epsilon_He3[7],color='y',s=size7)#Arnold 1978
plt.scatter(XS_Data_Q2eff_He3[3],XS_Data_Epsilon_He3[3],color='c',s=size3)#Dunn 1983
plt.scatter(XS_Data_Q2eff_He3[0],XS_Data_Epsilon_He3[0],color='b',s=size0)#Amroun 1994
plt.scatter(XS_Data_Q2eff_He3[5],XS_Data_Epsilon_He3[5],color='g',s=size5)#Nakagawa 2001
plt.scatter(XS_Data_Q2eff_He3[4],XS_Data_Epsilon_He3[4],color='m',s=size4)#Camsonne 2016
plt.scatter(XS_Data_Q2eff_He3[6],XS_Data_Epsilon_He3[6],color='darkorange',s=size6)#Barcus 2019
#plt.ylim([-1,1])

ax2 = ax.twiny()
ax2.scatter([], [],edgecolors='black',facecolors='none',s=outer_rad,label='0% Magnetic')
ax2.scatter([], [],color='black',edgecolors='black',s=outer_rad,label='100% Magnetic')

ax.legend(loc='best',fontsize=16,fancybox=True,framealpha=0.3)
ax2.legend(loc=(0.71,0.18),fontsize=16,fancybox=True,framealpha=0.3)
plt.show()


#Plot the residuals for each data point.
fig, ax = plt.subplots(figsize=(12,6))
ax.set_ylabel('Residual of SOG Fit',fontsize=18)
ax.set_xlabel(r'$Q^2$ (fm$^{-2}$)',fontsize=18)
ax.set_title(r'Residuals for the $^3$He Cross Section SOG Fit',fontsize=20)

#print(XS_Data_Q2eff_He3[0],XS_Data_ch_frac_He3[0])
plt.plot(XS_Data_Q2eff_He3[1],XS_Data_Residual_He3[1],'ro',label='Collard 1965')#Collard 1965
plt.plot(XS_Data_Q2eff_He3[2],XS_Data_Residual_He3[2],'o',color='saddlebrown',label='Szalata 1977')#Szalata 1977
plt.plot(XS_Data_Q2eff_He3[7],XS_Data_Residual_He3[7],'yo',label='Arnold 1978')#Arnold 1978
plt.plot(XS_Data_Q2eff_He3[3],XS_Data_Residual_He3[3],'co',label='Dunn 1983')#Dunn 1983
plt.plot(XS_Data_Q2eff_He3[0],XS_Data_Residual_He3[0],'bo',label='Amroun 1994')#Amroun 1994
plt.plot(XS_Data_Q2eff_He3[5],XS_Data_Residual_He3[5],'o',color='g',label='Nakagawa 2001')#Nakagawa 2001
plt.plot(XS_Data_Q2eff_He3[4],XS_Data_Residual_He3[4],'mo',label='Camsonne 2016')#Camsonne 2016
plt.plot(XS_Data_Q2eff_He3[6],XS_Data_Residual_He3[6],'o',color='darkorange',label='Barcus 2019')#Barcus 2019

#plt.ylim([-1,1])

ax.legend(loc='best',fontsize=16,fancybox=True,framealpha=0.3)
plt.show()

###############################################################################################
#Plot the charge and magnetic fractional contributions to the total XS for the 3H world data.#
###############################################################################################
ngaus = 8      #Number of Gaussians used to fit 3H data.
Z = 1          #Atomic number H3.

#My 3H thesis values.
R = (0.3,0.8,1.4,1.9,2.5,3.3,4.1,4.8)
Qich = (0.151488,0.348372,0.29635,0.0978631,0.121983,0.0242654,0.049329,4.40751e-11)
Qim = (0.190646,0.301416,0.318972,0.159433,0.173933,0.106361,0.0665564,0.0148866)

#Read in data.
with open('/home/skbarcus/JLab/SOG/3H_Data_Thesis_Datasets.txt') as f:
#with open('/home/skbarcus/JLab/SOG/3H_Data_Thesis.txt') as f:
    lines = f.readlines()

#Remove lines with column labels.
del lines[0]
del lines[0]

#Create arrays.
Raw_Data_H3 = []          #Array to store raw data.

#Read each line and split by the spaces and store as array.
#['Energy (GeV)', 'Theta (Degrees)', 'Sigma Experimental', 'Uncertainties', 'Dataset']
for line in lines:
    event = line.split()
    Raw_Data_H3.append(event)

#Turn data into numpy array and swap char to float.
Raw_Data_H3 = np.array(Raw_Data_H3)
Raw_Data_H3 = np.array(Raw_Data_H3.astype('float'))

#Examine data shape and check output.
#print('Raw_Data_H3.shape = ',Raw_Data_H3.shape)
#print('Raw_Data_H3[0] = ',Raw_Data_H3[0])

XS_Data_Q2eff_H3 = [[],[],[],[],[],[],[],[]]
XS_Data_ch_frac_H3 = [[],[],[],[],[],[],[],[]]
XS_Data_mag_frac_H3 = [[],[],[],[],[],[],[],[]]
XS_Data_Residual_H3 = [[],[],[],[],[],[],[],[]]
XS_Data_Epsilon_H3 = [[],[],[],[],[],[],[],[]]

for i in range(0,len(Raw_Data_H3)):
    #Calculate the residual per data point using the SOG fit result.
    res = (Raw_Data_H3[i][2] - XS_Data_H3(Raw_Data_H3[i][0],Raw_Data_H3[i][1])[0] ) /  XS_Data_H3(Raw_Data_H3[i][0],Raw_Data_H3[i][1])[0]
    #print('XSexp =',Raw_Data_He3[i][2],'   XSfit =',XS_Data_H3(Raw_Data_He3[i][0],Raw_Data_He3[i][1])[0],'   Residual =',res)
    if Raw_Data_H3[i][4]==1:#Collard 1965
        XS_Data_Q2eff_H3[0].append(XS_Data_H3(Raw_Data_H3[i][0],Raw_Data_H3[i][1])[5])
        XS_Data_ch_frac_H3[0].append(XS_Data_H3(Raw_Data_H3[i][0],Raw_Data_H3[i][1])[3])
        XS_Data_mag_frac_H3[0].append(XS_Data_H3(Raw_Data_H3[i][0],Raw_Data_H3[i][1])[4])
        XS_Data_Residual_H3[0].append(res)
        XS_Data_Epsilon_H3[0].append(XS_Data_H3(Raw_Data_H3[i][0],Raw_Data_H3[i][1])[8])
    if Raw_Data_H3[i][4]==2:#Beck 1984
        XS_Data_Q2eff_H3[1].append(XS_Data_H3(Raw_Data_H3[i][0],Raw_Data_H3[i][1])[5])
        XS_Data_ch_frac_H3[1].append(XS_Data_H3(Raw_Data_H3[i][0],Raw_Data_H3[i][1])[3])
        XS_Data_mag_frac_H3[1].append(XS_Data_H3(Raw_Data_H3[i][0],Raw_Data_H3[i][1])[4])
        XS_Data_Residual_H3[1].append(res)
        XS_Data_Epsilon_H3[1].append(XS_Data_H3(Raw_Data_H3[i][0],Raw_Data_H3[i][1])[8])
    if Raw_Data_H3[i][4]==3:#Amroun 1994
        XS_Data_Q2eff_H3[2].append(XS_Data_H3(Raw_Data_H3[i][0],Raw_Data_H3[i][1])[5])
        XS_Data_ch_frac_H3[2].append(XS_Data_H3(Raw_Data_H3[i][0],Raw_Data_H3[i][1])[3])
        XS_Data_mag_frac_H3[2].append(XS_Data_H3(Raw_Data_H3[i][0],Raw_Data_H3[i][1])[4])
        XS_Data_Residual_H3[2].append(res)
        XS_Data_Epsilon_H3[2].append(XS_Data_H3(Raw_Data_H3[i][0],Raw_Data_H3[i][1])[8])

XS_Data_Q2eff_H3 = np.array(XS_Data_Q2eff_H3)
XS_Data_ch_frac_H3= np.array(XS_Data_ch_frac_H3)
XS_Data_mag_frac_H3= np.array(XS_Data_mag_frac_H3)

#print(XS_Data_H3(Raw_Data_H3[0],Raw_Data_H3[1])[5])
#print(XS_Data_Q2eff_H3)
print('Raw_Data_H3.shape =',Raw_Data_H3.shape)
print('XS_Data_Q2eff_H3.shape =',XS_Data_Q2eff_H3.shape)

#Plot the world data by charge and magnetic XS fractions.
#Plot the charge contributions.
fig, ax = plt.subplots(figsize=(12,6))
ax.set_ylabel('Fractional Charge Contribution to Cross Section',fontsize=18)
ax.set_xlabel(r'$Q^2$ (fm$^{-2}$)',fontsize=18)
ax.set_title(r'Charge Contribution to the $^3$H Cross Section',fontsize=20)

#print(XS_Data_Q2eff_H3[0],XS_Data_ch_frac_H3[0])
plt.plot(XS_Data_Q2eff_H3[0],XS_Data_ch_frac_H3[0],'ro',label='Collard 1965')#Collard 1965
plt.plot(XS_Data_Q2eff_H3[1],XS_Data_ch_frac_H3[1],'go',label='Beck 1984')#Beck 1984
plt.plot(XS_Data_Q2eff_H3[2],XS_Data_ch_frac_H3[2],'bo',label='Amroun 1994')#Amroun 1994

ax.legend(loc='best',fontsize=16,fancybox=True,framealpha=0.3)
plt.show()

#Plot the magnetic contributions.
fig, ax = plt.subplots(figsize=(12,6))
ax.set_ylabel('Fractional Magnetic Contribution to Cross Section',fontsize=18)
ax.set_xlabel(r'$Q^2$ (fm$^{-2}$)',fontsize=18)
ax.set_title(r'Magnetic Contribution to the $^3$H Cross Section',fontsize=20)

#print(XS_Data_Q2eff_H3[0],XS_Data_ch_frac_H3[0])
plt.plot(XS_Data_Q2eff_H3[0],XS_Data_mag_frac_H3[0],'ro',label='Collard 1965')#Collard 1965
plt.plot(XS_Data_Q2eff_H3[1],XS_Data_mag_frac_H3[1],'go',label='Beck 1984')#Beck 1984
plt.plot(XS_Data_Q2eff_H3[2],XS_Data_mag_frac_H3[2],'bo',label='Amroun 1994')#Amroun 1994

ax.legend(loc='best',fontsize=16,fancybox=True,framealpha=0.3)
plt.show()

#Plot Q2eff vs Epsilon.
fig, ax = plt.subplots(figsize=(12,6))
ax.set_ylabel('Epsilon',fontsize=18)
ax.set_xlabel(r'$Q^2$ (fm$^{-2}$)',fontsize=18)
ax.set_title(r'World Data $Q^2$ vs. Epsilon',fontsize=20)

#print(XS_Data_Q2eff_H3[0],XS_Data_ch_frac_H3[0])
plt.plot(XS_Data_Q2eff_H3[0],XS_Data_Epsilon_H3[0],'ro',label='Collard 1965')#Collard 1965
plt.plot(XS_Data_Q2eff_H3[1],XS_Data_Epsilon_H3[1],'go',label='Beck 1984')#Beck 1984
plt.plot(XS_Data_Q2eff_H3[2],XS_Data_Epsilon_H3[2],'bo',label='Amroun 1994')#Amroun 1994

ax.legend(loc='best',fontsize=16,fancybox=True,framealpha=0.3)
plt.show()

#Plot the residuals for each data point.
fig, ax = plt.subplots(figsize=(12,6))
ax.set_ylabel('Residual of SOG Fit',fontsize=18)
ax.set_xlabel(r'$Q^2$ (fm$^{-2}$)',fontsize=18)
ax.set_title(r'Residuals for the $^3$H Cross Section SOG Fit',fontsize=20)

plt.plot(XS_Data_Q2eff_H3[0],XS_Data_Residual_H3[0],'ro',label='Collard 1965')#Collard 1965
plt.plot(XS_Data_Q2eff_H3[1],XS_Data_Residual_H3[1],'go',label='Beck 1984')#Beck 1984
plt.plot(XS_Data_Q2eff_H3[2],XS_Data_Residual_H3[2],'bo',label='Amroun 1994')#Amroun 1994

#plt.ylim([-1,1])

ax.legend(loc='best',fontsize=16,fancybox=True,framealpha=0.3)
plt.show()
