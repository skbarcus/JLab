import numpy as np
import matplotlib.pyplot as plt

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
E0 = 1                    #3.356 Initial e- energy GeV.#2.216
Ef = 0. 

def XS(Q2eff):
    return Q2eff

#Define Q2eff range to plot.
minq2 = 1     #0.48 for E0=1GeV
maxq2 = 17    #10.075
nsteps = 10000     #100
Q2eff = np.linspace(minq2,maxq2,nsteps)

def Q2_2_theta(Q2eff):
    #return 2*np.arcsin(  np.power( (1/(4*np.power(E0,2.)*GeV2fm/Q2eff-2*E0/MtHe3)) , 0.5 )  ) * rad2deg #Converts correctly.
    print('Q2 =',Q2eff,' Theta =', 2*np.arcsin(  np.power( (1/(4*np.power(E0,2.)*GeV2fm/Q2eff-2*E0/MtHe3)) , 0.5 )  )  * rad2deg)
    return 2*np.arcsin(  np.power( (1/(4*np.power(E0,2.)*GeV2fm/Q2eff-2*E0/MtHe3)) , 0.5 )  )  * rad2deg

    #return Q2eff/4

def theta_2_Q2(theta):
    #return 4*E0*E0*np.power(np.sin(theta*deg2rad/2.0),2.0) * GeV2fm / (1.0+2.0*E0*np.power(np.sin(theta*deg2rad/2.0),2.0)/MtHe3) #Converts correctly.
    #print('theta =',theta,' Q2 =',4*E0*E0*np.power(np.sin(theta*deg2rad/2.0),2.0) * GeV2fm / (1.0+2.0*E0*np.power(np.sin(theta*deg2rad/2.0),2.0)/MtHe3))
    return 4*E0*E0*np.power(np.sin(theta*deg2rad/2.0),2.0) * GeV2fm / (1.0+2.0*E0*np.power(np.sin(theta*deg2rad/2.0),2.0)/MtHe3)
    #return theta
    #return theta #*5 -77


Q2eff_test = np.linspace(1,21,100)
"""
for x in Q2eff_test:
    #Q2_2_theta(x)
    theta_2_Q2(x)

print('*******************************')
"""

fig, ax1 = plt.subplots(figsize=(12,6))

ax1.plot(Q2eff, XS(Q2eff), color='red', alpha=1.)

ax2 = ax1.secondary_xaxis('top',functions=(Q2_2_theta,theta_2_Q2))
#ax2 = ax1.secondary_xaxis('top',functions=(theta_2_Q2,Q2_2_theta))

dstep = 2
plt.xticks(np.arange(minq2,maxq2,step=dstep))
#plt.xticks(np.arange(0,61.0,step=dstep))

plt.show()




"""
fig, ax1 = plt.subplots(figsize=(12,6))
Q2eff2 = np.linspace(0.1,60,1000)

ax1.plot(Q2eff2, Q2_2_theta(Q2eff2), color='red', alpha=1.)

plt.show()


fig, ax1 = plt.subplots(figsize=(12,6))
theta1 = np.linspace(0.1,180,1000)

ax1.plot(theta1, theta_2_Q2(theta1), color='red', alpha=1.)

plt.show()
"""
