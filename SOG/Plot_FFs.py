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
alpha = 0.0072973525664             #Fine structure constant.
C = 299792458.0                     #Speed of light [m/s].\n",
e = 1.60217662E-19                  #Electron charge [C].\n",
e2_nuclear = 1.4399643929e-3        #Electron charge squared in nuclear units [GeV * fm].
MtHe3 = 3.0160293*0.9315            #Mass of He3 in GeV.
muHe3 = -2.1275*(3.0/2.0)           #Magnetic moment of 3He.
Z = 2                               #Atomic number of 3He.
A = 3                               #Number of nucleons.
gamma = 0.8*np.power(2.0/3.0,0.5);  #Gaussian width [fm] from Amroun gamma*sqrt(3/2) = 0.8 fm.
ngaus = 12

#My 3He thesis values.
R_He3_thesis = (0.3, 0.7, 0.9, 1.1, 1.5, 1.6, 2.2, 2.7, 3.3, 4.2, 4.3, 4.8)
Qich_He3_thesis = (0.0996392,0.214304,0.0199385,0.195676,0.0785533,0.167223,0.126926,0.0549379,0.0401401,0.0100803,0.0007217,4.98962e-12)
Qim_He3_thesis = (0.159649,0.0316168,0.277843,0.0364955,0.0329718,0.233469,0.117059,0.0581085,0.0485212,1.77602e-12,0.0240927,8.94934e-12)

#My 3H thesis values.
R_H3_thesis = (0.3,0.8,1.4,1.9,2.5,3.3,4.1,4.8)
Qich_H3_thesis = (0.151488,0.348372,0.29635,0.0978631,0.121983,0.0242654,0.049329,4.40751e-11)
Qim_H3_thesis = (0.190646,0.301416,0.318972,0.159433,0.173933,0.106361,0.0665564,0.0148866)

#Amroun 3He fit values.
R_He3_Amroun = (0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.,4.6,5.2) #Amroun Fit
Qich_He3_Amroun = (0.027614,0.170847,0.219805,0.170486,0.134453,0.100953,0.074310,0.053970,0.023689,0.017502,0.002034,0.004338)
Qim_He3_Amroun = (0.059785,0.138368,0.281326,0.000037,0.289808,0.019056,0.114825,0.042296,0.028345,0.018312,0.007843,0.)

#Amroun 3H fit values.
R_H3_Amroun = (0.1,0.5,0.9,1.3,1.6,2.0,2.4,2.9,3.4,4.,4.6,5.2) #Amroun Fit
Qich_H3_Amroun = (0.054706, 0.172505, 0.313852, 0.072056, 0.225333, 0.020849, 0.097374, 0.022273, 0.011933, 0.009121)#Amroun 3H
Qim_H3_Amroun = (0.075234, 0.164700, 0.273033, 0.037591, 0.252089, 0.027036, 0.098445, 0.040160, 0.016696, 0.015077)#Amroun 3H

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

#Read in the 3He data line by line.
with open('/home/skbarcus/JLab/SOG/New_Fits/3He_Fm_Data.txt') as f:
    lines = f.readlines()

#Remove first line with column labels.
del lines[0]
del lines[0]

#Create array to store data.
FF_Data = []

#Read each line and split by the spaces and store as array.
for line in lines:
    event = line.split()
    FF_Data.append(event)
    #print(event)

#Turn data into numpy array.
FF_Data = np.array(FF_Data)
FF_Data = FF_Data.astype('float')

#Examine data shape and check output.
print('FF_Data.shape = ',FF_Data.shape)
print('FF_Data[0] = ',FF_Data[0])

#Plot the residuals for each data point.
fig, ax = plt.subplots(figsize=(12,6))
ax.set_ylabel('|$F_m$|',fontsize=18)
ax.set_xlabel(r'$Q^2$ (fm$^{-2}$)',fontsize=18)
ax.set_title(r'$^3$He Magnetic Form Factor',fontsize=20)
ax.set_yscale('log')
#plt.ylim([-1,1])

#Define Q2eff range to plot.
Q2eff = np.linspace(0.00001,60,600)

#Plot Amroun and Representative fits.
plt.plot(Q2eff, np.absolute(Fm(Q2eff,Qim_He3_Amroun,R_He3_Amroun)), color='blue',label='Amroun Fit')
plt.plot(Q2eff, np.absolute(Fm(Q2eff,Qim_He3_thesis,R_He3_thesis)), color='black',label='Representative Fit')

FF_Data_Dataset = [[],[],[],[],[],[],[],[],]

#Sort by datasets.
for i in range(0,len(FF_Data)):
    if FF_Data[i][3]==2:#Collard 1964
        FF_Data_Dataset[0].append(FF_Data[i])
    if FF_Data[i][3]==11:#Bernheim 1977
        FF_Data_Dataset[1].append(FF_Data[i])
    if FF_Data[i][3]==9:#McCarthy 1977
        FF_Data_Dataset[2].append(FF_Data[i])
    if FF_Data[i][3]==4:#Dunn 1983
        FF_Data_Dataset[3].append(FF_Data[i])
    if FF_Data[i][3]==10:#Otterman 1984
        FF_Data_Dataset[4].append(FF_Data[i])
    if FF_Data[i][3]==1:#Amroun 1994
        FF_Data_Dataset[5].append(FF_Data[i])
    if FF_Data[i][3]==6:#Nakagawa 2001
        FF_Data_Dataset[6].append(FF_Data[i])
    if FF_Data[i][3]==5:#Camsonne 2016
        FF_Data_Dataset[7].append(FF_Data[i])

#Turn data into numpy array.
FF_Data_Dataset[0] = np.array(FF_Data_Dataset[0])
FF_Data_Dataset[1] = np.array(FF_Data_Dataset[1])
FF_Data_Dataset[2] = np.array(FF_Data_Dataset[2])
FF_Data_Dataset[3] = np.array(FF_Data_Dataset[3])
FF_Data_Dataset[4] = np.array(FF_Data_Dataset[4])
FF_Data_Dataset[5] = np.array(FF_Data_Dataset[5])
FF_Data_Dataset[6] = np.array(FF_Data_Dataset[6])
FF_Data_Dataset[7] = np.array(FF_Data_Dataset[7])
FF_Data_Dataset = np.array(FF_Data_Dataset)

#Examine data shape and check output.
print('FF_Data_Dataset.shape =',FF_Data_Dataset.shape)
print('FF_Data_Dataset[0][:,0] =',FF_Data_Dataset[0][:,0])
print('FF_Data_Dataset[0][:,1] =',FF_Data_Dataset[0][:,1])
print('FF_Data_Dataset[0][:,2] =',FF_Data_Dataset[0][:,2])
#print('FF_Data[:,0] =',FF_Data[:,0])
#print('FF_Data_Dataset = ',FF_Data_Dataset)

#Plot the world data.
#plt.errorbar(FF_Data[:,0],FF_Data[:,1],yerr = FF_Data[:,2],fmt='o',color='black',label='World Data')
plt.errorbar(FF_Data_Dataset[0][:,0],FF_Data_Dataset[0][:,1],yerr = FF_Data_Dataset[0][:,2],fmt='o',ms=5,color='red',label='Collard 1965')
plt.errorbar(FF_Data_Dataset[1][:,0],FF_Data_Dataset[1][:,1],yerr = FF_Data_Dataset[1][:,2],fmt='o',ms=5,color='saddlebrown',label='Bernheim 1977')
plt.errorbar(FF_Data_Dataset[2][:,0],FF_Data_Dataset[2][:,1],yerr = FF_Data_Dataset[2][:,2],fmt='o',ms=5,color='y',label='McCarthy 1977')
plt.errorbar([],[],0,label=' ',alpha=0)#Fake data to make blank space to arrange legend nicely around the data.
plt.errorbar(FF_Data_Dataset[3][:,0],FF_Data_Dataset[3][:,1],yerr = FF_Data_Dataset[3][:,2],fmt='o',ms=5,color='c',label='Dunn 1983')
plt.errorbar(FF_Data_Dataset[4][:,0],FF_Data_Dataset[4][:,1],yerr = FF_Data_Dataset[4][:,2],fmt='o',ms=5,color='purple',label='Otterman 1984')
plt.errorbar(FF_Data_Dataset[5][:,0],FF_Data_Dataset[5][:,1],yerr = FF_Data_Dataset[5][:,2],fmt='o',ms=5,color='blue',label='Amroun 1994')
plt.errorbar(FF_Data_Dataset[6][:,0],FF_Data_Dataset[6][:,1],yerr = FF_Data_Dataset[6][:,2],fmt='o',ms=5,color='g',label='Nakagawa 2001')
plt.errorbar(FF_Data_Dataset[7][:,0],FF_Data_Dataset[7][:,1],yerr = FF_Data_Dataset[7][:,2],fmt='o',ms=5,color='m',label='Camsonne 2016')
#Plot just my data point. Add errors of the ensemble 1-sigma band later.
plt.errorbar(34.19,2.968E-4,0,fmt='o',ms=5,color='darkorange',label='Barcus 2019')

ax.legend(loc='best',fontsize=16,fancybox=True,framealpha=0.3,ncol=2)
plt.show()
