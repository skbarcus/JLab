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

show_theory = 0
show_amroun = 0
show_rep = 1
show_ensemble = 1

He3_x2_cut = 437.312    #500 removes nonphysical fits (thesis). 437.312 = 1 sigma (582 out of 852 fits).
H3_x2_cut = 602.045     #603 removes nonphysical fits (thesis). 602.045 = 1 sigma (620 out of 908 fits).

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

#Read in the 3He data line by line.
with open('/home/skbarcus/JLab/SOG/New_Fits/Fch_Fit_Pars.txt') as f:
    lines = f.readlines()

#Remove first line with column labels.
del lines[0]

#Create array to store data.
Fch_He3_Fits = []

#Read each line and split by the spaces and store as array.
for line in lines:
    event = line.split()
    Fch_He3_Fits.append(event)
    #print(event)

#Turn data into numpy array.
Fch_He3_Fits = np.array(Fch_He3_Fits)

#Order of 3He fits array entries 0-41. R[0]=6, Q0ch=18 , Q0m=30.
#Chi2   rChi2   BIC   AIC    Qichtot   R[0]  R[1]  R[2]  R[3]  R[4]  R[5]  R[6]  R[7]  R[8]  R[9]  R[10]  R[11]  Q0ch    Q1ch    Q2ch    Q3ch    Q4ch    Q5ch    Q6ch    Q7ch    Q8ch    Q9ch    Q10ch    Q11ch

#Examine data shape and check output.
print('Fch_He3_Fits =',Fch_He3_Fits)
print('Fch_He3_Fits.shape = ',Fch_He3_Fits.shape)
print('Fch_He3_Fits[0] = ',Fch_He3_Fits[0])
