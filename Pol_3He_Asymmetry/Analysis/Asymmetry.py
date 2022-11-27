from matplotlib import pyplot as plt
import numpy as np
import math
import time
from matplotlib.transforms import Transform
from matplotlib.ticker import(AutoLocator, AutoMinorLocator)
from collections import OrderedDict

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

E0 = 2.18 #2.08 #Beam energy GeV.
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
    #print("Angles (degrees) = ",x)
    #print("Asymm Results = ",asymm)
    combined = OrderedDict(zip(x,asymm))
    #combined = set(combined)
    print("(Angles, Asymmetries) = ",combined)
    return np.array(asymm)

axtheta = 30
x = np.linspace(1,axtheta,500)
#x = np.linspace(18.3,18.8,50)
#x = [11.7]
#x = [11,13,15,16.3,17.5,18.5]
#x = [15, 17.13, 18.41, 19.45]
#x = [15, 17.13, 17.77,18.41]
#x = [8.5,11,13,13.5,14.5,16,16.4,18,18.5,20,22,24]
#x = [8.5,11,13,14.5,16,18,18.5,20,22,24,8.5,11,13.5,16.4,18.5,20,22,24]

axtheta = axtheta*pi/180
Q2max = 4*E0*Ef*pow(np.sin(axtheta/2),2)*GeV2fm
Q2effmax = Q2max*pow(1+(1.5*Z*alpha)/(E0*1.12*pow(A,1./3.)),2)
yax2 = []
for i in range(0,len(x)):
    yax2.append(i)

y = asymmetry(E0,x)

"""
shms_angles = (11,13,15,17.13,19.45)#High kinematic point.
shms_asymms = (0.0440913409471518,0.0644137073932522,0.0954473518134784,0.138384487219908,-0.154179671220445)#High kinematic point.
shms_asymm_err = (0.0039,0.0051,0.0173,0.0175,0.0331)#High kinematic point. #Didn't scale elastics by 4.7 for collimator (0.0018,0.0023,0.0079,0.0081,0.0153)
shms_times = ("1 Hour\n $11^{\circ}$","1 Hour\n $13^{\circ}$","1 Hour\n $15^{\circ}$","5 Hours\n $17.13^{\circ}$","16 Hours\n $19.45^{\circ}$")#High kinematic point.
"""
"""
shms_angles = (11,13,15,17.13,17.77)#New lower kinematic point.
shms_asymms = (0.0440913409471518,0.0644137073932522,0.0954473518134784,0.138384487219908,0.11392821942502)#New lower kinematic point.
shms_asymm_err = (0.0039,0.0051,0.0122,0.0138,0.0163)#New lower kinematic point.
shms_times = ("1 Hour\n $11^{\circ}$","1 Hour\n $13^{\circ}$","2 Hour\n $15^{\circ}$","8 Hours\n $17.13^{\circ}$","12 Hours\n $17.77^{\circ}$")#New lower kinematic point.
shms_angles = (11,13,15,17.13,17.77)#New lower kinematic point.
"""
"""
shms_angles = (11,13,15,17.1)#Hopefully final points.
shms_asymms = (0.0440913409471518,0.0644137073932522,0.0954473518134784,0.13835458900231506)#Hopefully final points.
shms_asymm_err = (0.0039,0.0062,0.0071,0.0071)#Hopefully final points.
shms_times = ("2 Hours\n $11^{\circ}$","2 Hours\n $13^{\circ}$","4 Hours\n $15^{\circ}$","16 Hours\n $17.1^{\circ}$")#Hopefully final points.
"""
shms_angles = (8.5,13)#Hopefully really final points.
shms_asymms = (0.026944154041414273,0.0692876293520045)#Hopefully really final points.
shms_asymm_err = (0.00072757,0.00253372)#Hopefully really final points.
shms_times = ("      8.5$^{\circ}$\n ~11 Hours Beam","      13$^{\circ}$\n ~6 Hours Beam")#Hopefully really final points.
hms_angles = (11.7,17)
hms_asymms = (0.053568445522965434,0.10381899441180453)
hms_asymm_err = (0.001317997,0.018068424)
hms_times = ("      11.7$^{\circ}$\n ~11 Hours Beam","      17$^{\circ}$\n ~6.5 Hours Beam")

#Add preliminary results. So far only SHMS.
shms_angles_prelim = (8.5,13)
shms_asymms_prelim = (0.02,0.0445)
shms_asymm_err_prelim = (0.0027,0.0041)
shms_times_prelim = ("      8.5$^{\circ}$\n ~11 Hours Beam","      13$^{\circ}$\n ~6 Hours Beam")#Hopefully really final points.
hms_angles_prelim = (11.7,17)
hms_asymms_prelim = (0.053568445522965434,0.10381899441180453)
hms_asymm_err_prelim = (0.001317997,0.018068424)
hms_times_prelim = ("      11.7$^{\circ}$\n ~11 Hours Beam","      17$^{\circ}$\n ~6.5 Hours Beam")
"""
hms_angles.append(18.41)
hms_asymms.append(-0.000458690566384221)
hms_asymm_err.append(0.0182)#Didn't scale elastics by 4.7 for collimator hms_asymm_err.append(0.0084)
hms_times.append("24 Hours\n $18.41^{\circ}$")
"""
"""
hms_angles.append(18.4)#Hopefully final points.
hms_asymms.append(0.0020174394461255816)#Hopefully final points.
hms_asymm_err.append(0.0202)#Hopefully final points.
hms_times.append("24 Hours\n $18.4^{\circ}$")#Hopefully final points.
"""
#hms_angles.append(17)#Hopefully final points.
#hms_asymms.append(0.10381899441180453)#Hopefully final points.
#hms_asymm_err.append(0.0212)#Hopefully final points.
#hms_times.append("24 Hours\n $17.5^{\circ}$")#Hopefully final points.

fig, ax = plt.subplots(figsize=(10,10))
plt.plot(x,y, 'red',linewidth=3,label="Predicted Asymmetry")#ax.plot(x,y, 'g') #Also works here.
plt.title("$^3$He Asymmetry for %.2f GeV" % E0,fontsize=36)
plt.xlabel("Scattering Angle (degrees)",fontsize=28)
plt.xticks(fontsize=24)
plt.ylabel("Physics Asymmetry",fontsize=28)
plt.yticks(fontsize=24)
ax.xaxis.set_major_locator(plt.MaxNLocator(15))

plt.errorbar(shms_angles_prelim, shms_asymms_prelim, yerr=shms_asymm_err_prelim,fmt='o',color='orange',label='SHMS Preliminary Asymmertries',linewidth=4,capsize=4,capthick=4,markersize=8)
ax.annotate(shms_times_prelim[0], (shms_angles_prelim[0], shms_asymms_prelim[0]), xytext=(shms_angles_prelim[0]-4.3, shms_asymms_prelim[0]-0.015),fontsize=18)
ax.annotate(shms_times_prelim[1], (shms_angles_prelim[1], shms_asymms_prelim[1]), xytext=(shms_angles_prelim[1]-1.0, shms_asymms_prelim[1]-0.017),fontsize=18)
#ax.annotate(shms_times_prelim[0], (shms_angles_prelim[0], shms_asymms_prelim[0]), xytext=(shms_angles_prelim[0]-3, shms_asymms_prelim[0]-0.013),fontsize=18)
#ax.annotate(shms_times_prelim[1], (shms_angles_prelim[1], shms_asymms_prelim[1]), xytext=(shms_angles_prelim[1]-2.5, shms_asymms_prelim[1]+0.007),fontsize=18)

#plt.errorbar(shms_angles, shms_asymms, yerr=shms_asymm_err,fmt='o',color='b',label='SHMS',linewidth=4,capsize=4,capthick=4,markersize=8)
plt.errorbar(hms_angles, hms_asymms, yerr=hms_asymm_err,fmt='o',color='g',label='HMS Projected (Rough Stat. Uncertainty)',linewidth=4,capsize=4,capthick=4,markersize=8)
#ax.annotate(shms_times[0], (shms_angles[0], shms_asymms[0]), xytext=(shms_angles[0]-3, shms_asymms[0]-0.013),fontsize=18)
#ax.annotate(shms_times[1], (shms_angles[1], shms_asymms[1]), xytext=(shms_angles[1]-2.5, shms_asymms[1]+0.007),fontsize=18)
ax.annotate(hms_times[0], (hms_angles[0], hms_asymms[0]), xytext=(hms_angles[0]-4, hms_asymms[0]+0.007),fontsize=18)
ax.annotate(hms_times[1], (hms_angles[1], hms_asymms[1]), xytext=(hms_angles[1]-0, hms_asymms[1]+0.0035),fontsize=18)

#Old label positions.
#ax.annotate(shms_times[0], (shms_angles[0], shms_asymms[0]), xytext=(shms_angles[0]-3, shms_asymms[0]-0.013),fontsize=18)
#ax.annotate(shms_times[1], (shms_angles[1], shms_asymms[1]), xytext=(shms_angles[1]-2.5, shms_asymms[1]+0.007),fontsize=18)
#ax.annotate(hms_times[0], (hms_angles[0], hms_asymms[0]), xytext=(hms_angles[0]-3, hms_asymms[0]-0.013),fontsize=18)
#ax.annotate(hms_times[1], (hms_angles[1], hms_asymms[1]), xytext=(hms_angles[1]-0, hms_asymms[1]+0.0035),fontsize=18)

#for i in range(len(shms_times)):
#    ax.annotate(shms_times[i], (shms_angles[i], shms_asymms[i]), xytext=(shms_angles[i]-1.75, shms_asymms[i]+0.0035),fontsize=20)
#for i in range(len(hms_times)):
#    ax.annotate(hms_times[i], (hms_angles[i], hms_asymms[i]), xytext=(hms_angles[i]-1.75, hms_asymms[i]+0.0035),fontsize=20)
#ax.annotate(shms_times, (shms_angles, shms_asymms), xytext=(shms_angles+0.05, shms_asymms+0.3), arrowprops=dict(facecolor='red', shrink=0.05))
plt.legend(loc="upper left",fontsize=15)
#ax1 = ax.twiny()
#ax1.set_xlabel('Q^2 fm^-2')
#ax1.set_xlim(0, Q2effmax)
plt.show()

#secax = ax.secondary_xaxis('top',functions=(x*100,x*0.01))#, functions=(x, yax2))
#secax.set_xlabel('angle [rad]')

print(type(x))
x.clear()
x.append(18.408)  #Min ~18.408 degrees.
y = asymmetry(E0,x)
print(y)
print("Asymmetry at %f degrees = %f " % (x[0],y[0]))

print("The script took %.2f seconds (%.2f minutes or %.2f hours) to run." % (time.time() - start, (time.time() - start)/60.,(time.time() - start)/60./60.)) #Print time to run.

