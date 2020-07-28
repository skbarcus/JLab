from matplotlib import pyplot as plt
import numpy as np
import math
import time
from matplotlib.transforms import Transform
from matplotlib.ticker import(AutoLocator, AutoMinorLocator)

start = time.time() #Start a timer.

pi = np.pi
deg2rad = pi/180.0 
GeV2fm = 1./0.0389              #Convert Q^2 units from GeV^2 to fm^-2.
hbar = 6.582e-16    #hbar in [eV*s].
C = 299792458.0                 #Speed of light [m/s]. 
e_ch = 1.60217662e-19              #Electron charge C.
alpha = 0.0072973525664         #1.0/137.0;              //Fine structure constant.
A = 3
Z = 2
muHe3 = -2.1275*(3.0/2.0)       #Dien's has this 3/2 factor for some reason, but it fits the data much better.  #2*2.793-1.913 is too naive.
MtHe3 = 3.0160293*0.9315        #Mass of He3 in GeV.
gamma = 0.8*pow(2.0/3.0,0.5)    #Gaussian width [fm] from Amroun gamma*sqrt(3/2) = 0.8 fm.
theta_pol = pi/2                #Polar polarization vector of target.
phi_pol = 0                     #Azimuthal polarization vector of target.

E0 = 2.08 #Beam energy GeV.
theta0 = 12 #Scattering angle in degrees.
Ef = E0/(1.0+2.0*E0*pow(np.sin(theta0/2.0),2.0)/MtHe3)
Q2 = 4*E0*Ef*pow(np.sin(theta0/2),2)*GeV2fm               #Converts energy in GeV^2 to fm^-2
Q2eff = Q2*pow(1+(1.5*Z*alpha)/(E0*1.12*pow(A,1./3.)),2)

#My 3He SOG FF fit parameters.
R = (0.3, 0.7, 0.9, 1.1, 1.5, 1.6, 2.2, 2.7, 3.3, 4.2, 4.3, 4.8)
Qich = (0.0996392,0.214304,0.0199385,0.195676,0.0785533,0.167223,0.126926,0.0549379,0.0401401,0.0100803,0.0007217,4.98962e-12)
Qim = (0.159649,0.0316168,0.277843,0.0364955,0.0329718,0.233469,0.117059,0.0581085,0.0485212,1.77602e-12,0.0240927,8.94934e-12)


def Fch(E,theta):
    #theta = theta*deg2rad
    #theta = x
    E0 = E
    Ef = E0/(1.0+2.0*E0*pow(np.sin(theta/2.0),2.0)/MtHe3)
    Q2 = 4*E0*Ef*pow(np.sin(theta/2),2)*GeV2fm               #Converts energy in GeV^2 to fm^-2
    Q2eff = Q2*pow(1+(1.5*Z*alpha)/(E0*1.12*pow(A,1./3.)),2)
    fitch = 0

    for i in range(0,len(R)):
        sumchtemp = (Qich[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( np.cos(pow(Q2eff,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (np.sin(pow(Q2eff,0.5)*R[i])/(pow(Q2eff,0.5)*R[i])) )
        fitch = fitch + sumchtemp
    fitch = fitch * np.exp(-0.25*Q2eff*pow(gamma,2.0))
    #print("Theta = %.2f degrees = %.3f radians, E0 = %.3f, Ef = %.3f, Q2 = %.3f, Q2eff = %.3f, fitch = %f " % (theta/deg2rad,theta,E0,Ef,Q2,Q2eff,fitch))
    return fitch

def Fm(E,theta):
    #theta = theta*deg2rad
    #theta = x
    E0 = E
    Ef = E0/(1.0+2.0*E0*pow(np.sin(theta/2.0),2.0)/MtHe3)
    Q2 = 4*E0*Ef*pow(np.sin(theta/2),2)*GeV2fm               #Converts energy in GeV^2 to fm^-2
    Q2eff = Q2*pow(1+(1.5*Z*alpha)/(E0*1.12*pow(A,1./3.)),2)
    fitm = 0

    for i in range(0,len(R)):
        summtemp = (Qim[i]/(1.0+2.0*pow(R[i],2.0)/pow(gamma,2.0))) * ( np.cos(pow(Q2eff,0.5)*R[i]) + (2.0*pow(R[i],2.0)/pow(gamma,2.0)) * (np.sin(pow(Q2eff,0.5)*R[i])/(pow(Q2eff,0.5)*R[i])) )
        fitm = fitm + summtemp
    fitm = fitm * np.exp(-0.25*Q2eff*pow(gamma,2.0))
   # print("Theta = %.2f degrees = %.3f radians, E0 = %.3f, Ef = %.3f, Q2 = %.3f, Q2eff = %.3f, fitm = %f " % (theta/deg2rad,theta,E0,Ef,Q2,Q2eff,fitm))
    return fitm

def asymmetry(E,x):
    asymm = [] #List to store asymmetry results

    for xval in x:
        theta = xval
        #print("xval = ",xval)
        theta = theta*deg2rad       #Convert degrees to radians.
        E0 = E
        Ef = E0/(1.0+2.0*E0*pow(np.sin(theta/2.0),2.0)/MtHe3)
        Q2 = 4*E0*Ef*pow(np.sin(theta/2),2)*GeV2fm               #Converts energy in GeV^2 to fm^-2
        Q2eff = Q2*pow(1+(1.5*Z*alpha)/(E0*1.12*pow(A,1./3.)),2)
        tau = Q2eff/(4*pow(MtHe3,2.)*GeV2fm)  #Tau using fm^-2.
        epsilon = pow( 1 + 2*(1+tau)*pow(np.tan(theta),2.) ,-1.)
        
        asymm.append( ( -2 * pow(tau * (1+tau),0.5) * np.tan(theta/2) ) / ( pow(Fch(E0,theta),2.) + (tau/epsilon) * (pow(muHe3,2.)) * pow(Fm(E0,theta),2.) ) * ( np.sin(theta_pol)*np.cos(phi_pol)*Fch(E0,theta)*Fm(E0,theta)*muHe3 + pow(tau*(1+(1+tau)*pow(np.tan(theta/2),2.)),0.5)*np.cos(theta_pol)*pow(Fm(E0,theta),2.)*pow(muHe3,2.) ) )
        #print("Theta = %.2f degrees = %.3f radians, E0 = %.3f, Ef = %.3f, Q2 = %.3f fm^-2, Q2eff = %.3f fm^-2, tau = %f, epsilon = %f " % (theta/deg2rad,theta,E0,Ef,Q2,Q2eff,tau,epsilon))
        #print("Asymm Results = ",asymm)
    return np.array(asymm)

axtheta = 30
x = np.linspace(1,axtheta,500)
axtheta = axtheta*pi/180
Q2max = 4*E0*Ef*pow(np.sin(axtheta/2),2)*GeV2fm
Q2effmax = Q2max*pow(1+(1.5*Z*alpha)/(E0*1.12*pow(A,1./3.)),2)
yax2 = []
for i in range(0,len(x)):
    yax2.append(i)

y = asymmetry(E0,x)

fig, ax = plt.subplots(figsize=(10,10))
plt.plot(x,y, 'red')#ax.plot(x,y, 'g') #Also works here.
plt.title("3He Asymmetry for %.2f GeV" % E0)
plt.xlabel("Scattering Angle (degrees)")
plt.ylabel("Asymmetry")
ax.xaxis.set_major_locator(plt.MaxNLocator(30))
#ax1 = ax.twiny()
#ax1.set_xlabel('Q^2 fm^-2')
#ax1.set_xlim(0, Q2effmax)
plt.show()

#secax = ax.secondary_xaxis('top',functions=(x*100,x*0.01))#, functions=(x, yax2))
#secax.set_xlabel('angle [rad]')

x = []
x.append(18.408)  #Min ~18.408 degrees.
y = asymmetry(E0,theta0)
print(y)
print("Asymmetry at %f degrees = %f " % (x[0],y[0]))

print("The script took %.2f seconds (%.2f minutes or %.2f hours) to run." % (time.time() - start, (time.time() - start)/60.,(time.time() - start)/60./60.)) #Print time to run.

