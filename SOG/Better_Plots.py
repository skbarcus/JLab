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

show_theory = 0
show_amroun = 1
show_rep = 1
show_ensemble = 1

He3_x2_cut = 500    #500 removes nonphysical fits (thesis). 437.312 = 1 sigma (582 out of 852 fits).
H3_x2_cut = 750     #603 removes nonphysical fits (thesis). 602.045 = 1 sigma (620 out of 908 fits).

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

#Create arrays to hold radii values.
radii_3He = []
radii_3H = []

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
#with open('/home/skbarcus/JLab/SOG/Ri_Fits_Final_n=12_1352_12_22_2018.txt') as f:
#with open('/home/skbarcus/JLab/SOG/Fits_3He_Sum1.txt') as f:
with open('/home/skbarcus/JLab/SOG/All_Fit_Pars_3He_4-13-2022.txt') as f:
    lines = f.readlines()

#Remove first line with column labels.
del lines[0]

#Create array to store data.
He3_Fits = []

#Read each line and split by the spaces and store as array.
for line in lines:
    event = line.split()
    He3_Fits.append(event)
    #print(event)

#Turn data into numpy array.
He3_Fits = np.array(He3_Fits)

#Order of 3He fits array entries 0-41. R[0]=6, Q0ch=18 , Q0m=30.
#Chi2   rChi2   BIC   AIC    Qichtot    Qimtot  R[0]  R[1]  R[2]  R[3]  R[4]  R[5]  R[6]  R[7]  R[8]  R[9]  R[10]  R[11]  Q0ch    Q1ch    Q2ch    Q3ch    Q4ch    Q5ch    Q6ch    Q7ch    Q8ch    Q9ch    Q10ch    Q11ch    Q0m    Q1m    Q2m    Q3m    Q4m    Q5m    Q6m    Q7m    Q8m    Q9m    Q10m    Q11m

#Examine data shape and check output.
print('He3_Fits.shape = ',He3_Fits.shape)
print('He3_Fits[0] = ',He3_Fits[0])

#for i in range(0,len(He3_Fits)):
    #print('Entry',i,' Length',len(He3_Fits[i]))

#Split array into targets (elastic or not elastic) and training data.
He3_Fits = np.hsplit(He3_Fits,np.array([6,18,30]))
print('type(He3_Fits[0])',type(He3_Fits[0]))
print('type(He3_Fits[0][0])',type(He3_Fits[0][0]))
#He3_Fits[0] = np.array(He3_Fits[0])
stats_He3 = np.array(He3_Fits[0].astype('float'))
R_He3 = np.array(He3_Fits[1].astype('float'))
Qich_He3 = np.array(He3_Fits[2].astype('float'))
Qim_He3 = np.array(He3_Fits[3].astype('float'))

print('stats_He3.shape',stats_He3.shape)
print('stats_He3[0]',stats_He3[0])

print('R_He3.shape',R_He3.shape)
print('R_He3[0]',R_He3[0])

print('Qich_He3.shape',Qich_He3.shape)
print('Qich_He3[0]',Qich_He3[0])

print('Qim_He3.shape',Qim_He3.shape)
print('Qim_He3[0]',Qim_He3[0])

#Read in the 3H data line by line. Remember last 4 entries for R, Qich, and Qim are meaningless and can just be ignored.
#with open('/home/skbarcus/JLab/SOG/Ri_Fits_3H_Final_n=8_2600_12_22_2018.txt') as f:
#with open('/home/skbarcus/JLab/SOG/Fits_3H_Sum1.txt') as f:
with open('/home/skbarcus/JLab/SOG/All_Fit_Pars_3H_4-13-2022.txt') as f:
    lines = f.readlines()

#Remove first line with column labels.
del lines[0]

#Create array to store data.
H3_Fits = []

#Read each line and split by the spaces and store as array.
for line in lines:
    event = line.split()
    H3_Fits.append(event)
    #print(event)

#Turn data into numpy array.
H3_Fits = np.array(H3_Fits)

#Order of 3He fits array entries 0-41. R[0]=6, Q0ch=18 , Q0m=30.
#Chi2   rChi2   BIC   AIC    Qichtot    Qimtot  R[0]  R[1]  R[2]  R[3]  R[4]  R[5]  R[6]  R[7]  R[8]  R[9]  R[10]  R[11]  Q0ch    Q1ch    Q2ch    Q3ch    Q4ch    Q5ch    Q6ch    Q7ch    Q8ch    Q9ch    Q10ch    Q11ch    Q0m    Q1m    Q2m    Q3m    Q4m    Q5m    Q6m    Q7m    Q8m    Q9m    Q10m    Q11m

#Examine data shape and check output.
print('H3_Fits.shape = ',H3_Fits.shape)
print('H3_Fits[0] = ',H3_Fits[0])

#Split array into targets (elastic or not elastic) and training data.
H3_Fits = np.hsplit(H3_Fits,np.array([6,14,18,26,30,38]))

stats_H3 = np.array(H3_Fits[0].astype('float'))
R_H3 = np.array(H3_Fits[1].astype('float'))
Qich_H3 = np.array(H3_Fits[3].astype('float'))
Qim_H3 = np.array(H3_Fits[5].astype('float'))

print('stats_H3.shape',stats_H3.shape)
print('stats_H3[0]',stats_H3[0])

print('R_H3.shape',R_H3.shape)
print('R_H3[0]',R_H3[0])

print('Qich_H3.shape',Qich_H3.shape)
print('Qich_H3[0]',Qich_H3[0])

print('Qim_H3.shape',Qim_H3.shape)
print('Qim_H3[0]',Qim_H3[0])

"""
#Read in theory curves.
with open('/home/skbarcus/Tritium/Analysis/SOG/3He_Fch_Conventional_Q2.txt') as f:
    lines = f.readlines()

#Create array to store data.
He3_fch_conv = []

print('type(lines)',type(lines))

#Read each line and split by the spaces and store as array.
for line in lines:
    #print('line',line)
    event = line.split()
    He3_fch_conv.append(event)
    #print(event)

#Turn data into numpy array.
He3_fch_conv = np.array(He3_fch_conv)

print('He3_fch_conv.shape',He3_fch_conv.shape)

#Split array into x and y.
He3_fch_conv = np.hsplit(He3_fch_conv,np.array([1]))

#Save x and y as numpy arrays of type float.
He3_fch_conv_x = np.array(He3_fch_conv[0].astype('float'))
He3_fch_conv_y = np.array(He3_fch_conv[1].astype('float'))

print('He3_fch_conv_x.shape',He3_fch_conv_x.shape)
print('He3_fch_conv_y.shape',He3_fch_conv_y.shape)
"""

#Read in all theory curves and build the x,y arrays for plotting.
files = ['/home/skbarcus/Tritium/Analysis/SOG/3He_Fch_Conventional_Q2.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3He_Fch_CST_Q2.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3He_Fch_XEFT500_Q2.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3He_Fch_XEFT600_Q2.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3H_Fch_Conventional_Q2.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3H_Fch_CST_Q2.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3H_Fch_XEFT500_Q2.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3H_Fch_XEFT600_Q2.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3He_Fm_Conventional_Q2.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3He_Fm_CST_Q2.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3He_Fm_XEFT500_Q2.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3He_Fm_XEFT600_Q2.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3H_Fm_Conventional_Q2.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3H_Fm_CST_Q2.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3H_Fm_XEFT500_Q2.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3H_Fm_XEFT600_Q2.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3He_Fch_Amroun_Error_Band_Down.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3He_Fch_Amroun_Error_Band_Up.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3He_Fm_Amroun_Error_Band_Down.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3He_Fm_Amroun_Error_Band_Up.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3H_Fch_Amroun_Error_Band_Down.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3H_Fch_Amroun_Error_Band_Up.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3H_Fm_Amroun_Error_Band_Down.txt',
         '/home/skbarcus/Tritium/Analysis/SOG/3H_Fm_Amroun_Error_Band_Up.txt']

lines = []

for file in range(0,len(files)):
    with open(files[file]) as f:
        #lines[file] = f.readlines()
        lines.append(f.readlines())

#print('lines',lines)
#lines = np.array(lines)
#print('lines.shape',lines.shape)

print('type(lines[0])',type(lines[0]))

#Create array to store data.
He3_fch_conv_test = []
He3_fch_CST_test = []
theory_test = []

#Read each line and split by the spaces and store as array.
for i in range(0,len(lines)):
    theory_test.append([])
    for j in range(0,len(lines[i])):
        line = lines[i][j]
        #print('line',line)
        event = line.split()
        theory_test[i].append(event)
        #theory_test.append(event)

#print('theory_test',theory_test)

#Turn data into numpy array.
theory_test = np.array(theory_test)
print('theory_test',theory_test.shape)

#Split array into x and y.
for i in range(0,len(theory_test)):
    theory_test[i] = np.array(theory_test[i])
    theory_test[i] = np.hsplit(theory_test[i],np.array([1]))

#print('theory_test',theory_test)
#print('theory_test[0][0]',theory_test[0][0])

#Save x and y as numpy arrays of type float.
He3_fch_conv_x = np.array(theory_test[0][0].astype('float'))
He3_fch_conv_y = np.array(theory_test[0][1].astype('float'))
He3_fch_CST_x = np.array(theory_test[1][0].astype('float'))
He3_fch_CST_y = np.array(theory_test[1][1].astype('float'))
He3_Fch_XEFT500_x = np.array(theory_test[2][0].astype('float'))
He3_Fch_XEFT500_y = np.array(theory_test[2][1].astype('float'))
He3_Fch_XEFT600_x = np.array(theory_test[3][0].astype('float'))
He3_Fch_XEFT600_y = np.array(theory_test[3][1].astype('float'))

H3_fch_conv_x = np.array(theory_test[4][0].astype('float'))
H3_fch_conv_y = np.array(theory_test[4][1].astype('float'))
H3_fch_CST_x = np.array(theory_test[5][0].astype('float'))
H3_fch_CST_y = np.array(theory_test[5][1].astype('float'))
H3_Fch_XEFT500_x = np.array(theory_test[6][0].astype('float'))
H3_Fch_XEFT500_y = np.array(theory_test[6][1].astype('float'))
H3_Fch_XEFT600_x = np.array(theory_test[7][0].astype('float'))
H3_Fch_XEFT600_y = np.array(theory_test[7][1].astype('float'))

He3_fm_conv_x = np.array(theory_test[8][0].astype('float'))
He3_fm_conv_y = np.array(theory_test[8][1].astype('float'))
He3_fm_CST_x = np.array(theory_test[9][0].astype('float'))
He3_fm_CST_y = np.array(theory_test[9][1].astype('float'))
He3_Fm_XEFT500_x = np.array(theory_test[10][0].astype('float'))
He3_Fm_XEFT500_y = np.array(theory_test[10][1].astype('float'))
He3_Fm_XEFT600_x = np.array(theory_test[11][0].astype('float'))
He3_Fm_XEFT600_y = np.array(theory_test[11][1].astype('float'))

H3_fm_conv_x = np.array(theory_test[12][0].astype('float'))
H3_fm_conv_y = np.array(theory_test[12][1].astype('float'))
H3_fm_CST_x = np.array(theory_test[13][0].astype('float'))
H3_fm_CST_y = np.array(theory_test[13][1].astype('float'))
H3_Fm_XEFT500_x = np.array(theory_test[14][0].astype('float'))
H3_Fm_XEFT500_y = np.array(theory_test[14][1].astype('float'))
H3_Fm_XEFT600_x = np.array(theory_test[15][0].astype('float'))
H3_Fm_XEFT600_y = np.array(theory_test[15][1].astype('float'))

He3_Fch_Amroun_Error_Band_Down_x = np.array(theory_test[16][0].astype('float'))
He3_Fch_Amroun_Error_Band_Down_y = np.array(theory_test[16][1].astype('float'))
He3_Fch_Amroun_Error_Band_Up_x = np.array(theory_test[17][0].astype('float'))
He3_Fch_Amroun_Error_Band_Up_y = np.array(theory_test[17][1].astype('float'))

He3_Fm_Amroun_Error_Band_Down_x = np.array(theory_test[18][0].astype('float'))
He3_Fm_Amroun_Error_Band_Down_y = np.array(theory_test[18][1].astype('float'))
He3_Fm_Amroun_Error_Band_Up_x = np.array(theory_test[19][0].astype('float'))
He3_Fm_Amroun_Error_Band_Up_y = np.array(theory_test[19][1].astype('float'))

H3_Fch_Amroun_Error_Band_Down_x = np.array(theory_test[20][0].astype('float'))
H3_Fch_Amroun_Error_Band_Down_y = np.array(theory_test[20][1].astype('float'))
H3_Fch_Amroun_Error_Band_Up_x = np.array(theory_test[21][0].astype('float'))
H3_Fch_Amroun_Error_Band_Up_y = np.array(theory_test[21][1].astype('float'))

H3_Fm_Amroun_Error_Band_Down_x = np.array(theory_test[22][0].astype('float'))
H3_Fm_Amroun_Error_Band_Down_y = np.array(theory_test[22][1].astype('float'))
H3_Fm_Amroun_Error_Band_Up_x = np.array(theory_test[23][0].astype('float'))
H3_Fm_Amroun_Error_Band_Up_y = np.array(theory_test[23][1].astype('float'))

print('He3_fch_conv_x.shape',He3_fch_conv_x.shape)
print('He3_fch_conv_y.shape',He3_fch_conv_y.shape)

def gaus_fit(x, C, mu, sigma):
    return ( C * np.exp(-1.0 * (x - mu)**2 / (2 * sigma**2)) )

#Define charge form factor function.
def Fch(Q2eff,Qich,R):
    sumFch_ff = 0
    Fch_ff = 0
    for i in range(ngaus):
        sumFch_ff = (Qich[i]/(1.0+2.0*np.power(R[i],2.0)/np.power(gamma,2.0))) * ( np.cos(np.power(Q2eff,0.5)*R[i]) + (2.0*np.power(R[i],2.0)/np.power(gamma,2.0)) * (np.sin(np.power(Q2eff,0.5)*R[i])/(np.power(Q2eff,0.5)*R[i])) )
        Fch_ff = Fch_ff + sumFch_ff
    Fch_ff =  Fch_ff * np.exp(-0.25*Q2eff*np.power(gamma,2.0))
    return Fch_ff

#Define charge form factor function for derivative calculation.
def Fch_deriv(Q2eff):
    sumFch_ff = 0
    Fch_ff = 0
    for i in range(ngaus):
        sumFch_ff = (Qich[i]/(1.0+2.0*np.power(Ri[i],2.0)/np.power(gamma,2.0))) * ( np.cos(np.power(Q2eff,0.5)*Ri[i]) + (2.0*np.power(Ri[i],2.0)/np.power(gamma,2.0)) * (np.sin(np.power(Q2eff,0.5)*Ri[i])/(np.power(Q2eff,0.5)*Ri[i])) )
        Fch_ff = Fch_ff + sumFch_ff
    Fch_ff =  Fch_ff * np.exp(-0.25*Q2eff*np.power(gamma,2.0))
    return Fch_ff

#Define function to calculate rms charge radius.
def rms_radius(deriv):
    radius = pow(-6*deriv,0.5)
    return radius

#Define magnetic form factor function.
def Fm(Q2eff,Qim,R):
    sumFm_ff = 0
    Fm_ff = 0
    for i in range(ngaus):
        sumFm_ff = (Qim[i]/(1.0+2.0*np.power(R[i],2.0)/np.power(gamma,2.0))) * ( np.cos(np.power(Q2eff,0.5)*R[i]) + (2.0*np.power(R[i],2.0)/np.power(gamma,2.0)) * (np.sin(np.power(Q2eff,0.5)*R[i])/(np.power(Q2eff,0.5)*R[i])) )
        Fm_ff = Fm_ff + sumFm_ff
    Fm_ff =  Fm_ff * np.exp(-0.25*Q2eff*np.power(gamma,2.0))
    return Fm_ff

#Define the charge density from I. Sick. 
def rho_ch(r,Qich,R,Z):
    rho = 0;
    rho_temp = 0;
   
    for i in range(ngaus):
        rho_temp = Qich[i]/( 1+2*np.power(R[i],2.)/np.power(gamma,2.) ) * (  np.exp( -np.power((r-R[i]),2.)/np.power(gamma,2.) ) + np.exp( -np.power((r+R[i]),2.)/np.power(gamma,2.) )  );
        rho = rho + rho_temp;
    
    rho = Z/(2*np.power(pi,1.5)*np.power(gamma,3.)) * rho; #Really Z*e factor but to make the units of rho be e/fm^3 I divided out e here.
    return rho;

#Define the charge density from I. Sick. 
def rho_ch_int(r,Qich,R,Z):
    rho = 0;
    rho_temp = 0;
   
    for i in range(ngaus):
        rho_temp = Qich[i]/( 1+2*np.power(R[i],2.)/np.power(gamma,2.) ) * (  np.exp( -np.power((r-R[i]),2.)/np.power(gamma,2.) ) + np.exp( -np.power((r+R[i]),2.)/np.power(gamma,2.) )  );
        rho = rho + rho_temp;
    
    rho = 4*pi*np.power(r,2.) * Z/(2*np.power(pi,1.5)*np.power(gamma,3.)) * rho; #Really Z*e factor but to make the units of rho be e/fm^3 I divided out e here.
    return rho;


x0 = 1

def f(q):
 return q*q+x0*q

d = derivative(f, 1.0, dx=1e-5)
print ('d=',d)

Qich = Qich_He3_thesis
Ri = R_He3_thesis


print('Fch_deriv(0.01)',Fch_deriv(0.013605))
d = derivative(Fch_deriv, 1.0, dx=1e-5)
print ('Fch_deriv=',d)

"""
#x = np.linspace(0.00001,60,600)
x = np.arange(0,4,0.01)
y = np.arange(0,4,0.01)

d = derivative(f(y), 1.0, dx=1e-3)
print ('d=',d)
"""

#Find 1-sigma (68.27%) band in the 3He ensemble of fits. 
#Create arrays to hold the sorted stat values to build uncertainty bands.
He3_x2 = []
He3_BIC = []
He3_AIC = []

#Fill the stat arrays.
for fit in range(0,len(stats_He3)):
    if stats_He3[fit][0]<He3_x2_cut:
        He3_x2.append(stats_He3[fit][0])
        He3_BIC.append(stats_He3[fit][2])
        He3_AIC.append(stats_He3[fit][3])

He3_x2 = np.array(He3_x2)
He3_BIC = np.array(He3_BIC)
He3_AIC = np.array(He3_AIC)

#Sort arrays from lowest value to highest value.
He3_x2 = np.sort(He3_x2)
He3_BIC = np.sort(He3_BIC)
He3_AIC = np.sort(He3_AIC)

#Calculate the 1-sigma cut values.
He3_nfits = len(He3_x2)
He3_nfits_1sig = round(len(He3_x2)*0.6827)
He3_x2_1sig = He3_x2[round(len(He3_x2)*0.6827)]
He3_BIC_1sig = He3_BIC[round(len(He3_BIC)*0.6827)]
He3_AIC_1sig = He3_AIC[round(len(He3_AIC)*0.6827)]

#Print 1-sigma results.
print('len(He3_x2) =',He3_nfits)
print('Number of He3 fits within 1-sigma band =',He3_nfits_1sig)
print('He3 1-sigma x2 =',He3_x2_1sig)
print('He3 1-sigma BIC =',He3_BIC_1sig)
print('He3 1-sigma AIC =',He3_AIC_1sig)

#Find 1-sigma (68.27%) band in the 3H ensemble of fits. 
#Create arrays to hold the sorted stat values to build uncertainty bands.
H3_x2 = []
H3_BIC = []
H3_AIC = []

#Fill the stat arrays.
for fit in range(0,len(stats_H3)):
    if stats_H3[fit][0]<H3_x2_cut:
        H3_x2.append(stats_H3[fit][0])
        H3_BIC.append(stats_H3[fit][2])
        H3_AIC.append(stats_H3[fit][3])

H3_x2 = np.array(H3_x2)
H3_BIC = np.array(H3_BIC)
H3_AIC = np.array(H3_AIC)

#Sort arrays from lowest value to highest value.
H3_x2 = np.sort(H3_x2)
H3_BIC = np.sort(H3_BIC)
H3_AIC = np.sort(H3_AIC)

#Calculate the 1-sigma cut values.
H3_nfits = len(H3_x2)
H3_nfits_1sig = round(len(H3_x2)*0.6827)
H3_x2_1sig = H3_x2[round(len(H3_x2)*0.6827)]
H3_BIC_1sig = H3_BIC[round(len(H3_BIC)*0.6827)]
H3_AIC_1sig = H3_AIC[round(len(H3_AIC)*0.6827)]

#Print 1-sigma results.
print('len(H3_x2) =',H3_nfits)
print('Number of H3 fits within 1-sigma band =',H3_nfits_1sig)
print('H3 1-sigma x2 =',H3_x2_1sig)
print('H3 1-sigma BIC =',H3_BIC_1sig)
print('H3 1-sigma AIC =',H3_AIC_1sig)

#Turn on TeX for labels.
plt.rcParams['text.usetex'] = True

#Plot the 3He charge FF Fits.
fig, ax = plt.subplots(figsize=(12,6))
ax.set_title('$^3$He Charge Form Factor',fontsize=20)
ax.set_ylabel('$F_{ch}(Q^2)$',fontsize=16)
ax.set_xlabel('$Q^2$ (fm$^{-2}$)',fontsize=16)
ax.set_yscale('log')

min = 0
max = 60
plt.xticks(np.arange(min, max+1, 2.0))

#Define Q2eff range to plot.
Q2eff = np.linspace(0.00001,60,600)

#derivative = Fch(Q2eff,Qich_He3_thesis,R_He3_thesis).deriv()

#d = derivative(Fch(Q2eff,Qich_He3_thesis,R_He3_thesis), 1.0, dx=1e-3)
#print ('d=',d)



#Plot ensemble of 3He fits surviving x^2 cut.
if show_ensemble==1:
    for fit in range(0,len(R_He3)):
        if stats_He3[fit][0]<He3_x2_cut:
            plt.plot(Q2eff, np.absolute(Fch(Q2eff,Qich_He3[fit],R_He3[fit])), color='red', alpha=0.2)
            Ri = R_He3[fit]
            Qich = Qich_He3[fit]
            d = derivative(Fch_deriv, 0.0015, dx=1e-5)
            radius = rms_radius(d)
            radii_3He.append(radius)
            #print ('Fch_deriv =',d,'radius =',radius)

radii_3He = np.array(radii_3He)

#Dummy plot just to add a label for the ensemble fits.
plt.plot(0, 0, color='red',label='Ensemble Fits')

#Plot Amroun charge FF representative fit.
if show_amroun==1:
    plt.plot(Q2eff, np.absolute(Fch(Q2eff,Qich_He3_Amroun,R_He3_Amroun)), color='blue',label='Amroun Fit')
    #This works to fill the error bands. Need to flip one error band array so the fill works correctly.
    x = np.append(He3_Fch_Amroun_Error_Band_Down_x,np.flip(He3_Fch_Amroun_Error_Band_Up_x))
    y = np.append(He3_Fch_Amroun_Error_Band_Down_y,np.flip(He3_Fch_Amroun_Error_Band_Up_y))
    plt.fill(x,y)

#Plot new 3He charge FF representative fit.
if show_rep==1:
    plt.plot(Q2eff, np.absolute(Fch(Q2eff,Qich_He3_thesis,R_He3_thesis)), color='black',label='New Representative Fit') #Plot 3He representative fit.

#Plot theory curves.
#Make curves smoother with spline. Doesn't help much. Leave for possible future use.
#He3_fch_conv_x = np.reshape(He3_fch_conv_x,(60,))
#He3_fch_conv_y = np.reshape(He3_fch_conv_y,(60,))
#X_Y_Spline = make_interp_spline(He3_fch_conv_x, He3_fch_conv_y)
#cubic_interploation_model = interp1d(He3_fch_conv_x, He3_fch_conv_y, kind = "cubic")
#X_ = Q2eff
#Y_ = X_Y_Spline(X_)
#Y_=cubic_interploation_model(X_)
#plt.plot(X_, Y_, color='green',label='Conventional Approach Marcucci 2016')

if show_theory==1:
    plt.plot(He3_fch_conv_x, He3_fch_conv_y, color='green',label='Conventional Approach Marcucci 2016')
    plt.plot(He3_fch_CST_x, He3_fch_CST_y, color='cyan',label='Covalent Spectator Theorem Marcucci 2016')
    plt.plot(He3_Fch_XEFT500_x, He3_Fch_XEFT500_y, color='m',label='$\chi$EFT500 Marcucci 2016')
    plt.plot(He3_Fch_XEFT600_x, He3_Fch_XEFT600_y, color='brown',label='$\chi$EFT600 Marcucci 2016')

#plt.plot(He3_Fch_Amroun_Error_Band_Down_x, He3_Fch_Amroun_Error_Band_Down_y, color='orange')
#plt.plot(He3_Fch_Amroun_Error_Band_Up_x, He3_Fch_Amroun_Error_Band_Up_y, color='orange')

"""Bad idea to try a poly fit on error bands.
He3_Fch_Amroun_Error_Band_Down_Fit = Polynomial.fit(He3_Fch_Amroun_Error_Band_Down_x, He3_Fch_Amroun_Error_Band_Down_y, deg=50)
plt.plot(*He3_Fch_Amroun_Error_Band_Down_Fit.linspace(),color='yellow')
"""

ax.legend(loc='upper right')
plt.show()

#Plot the 3He charge distributions.
fig, ax = plt.subplots(figsize=(12,6))
ax.set_title('$^3$He Charge Density',fontsize=20)
ax.set_ylabel(r'$\rho(r)$ (e/fm$^3$)',fontsize=16)
ax.set_xlabel('Radius (fm)',fontsize=16)
#ax.set_yscale('log')

#Define radii range to plot.
r = np.linspace(0.00001,5,5000)

if show_ensemble==1:
    for fit in range(0,len(R_He3)):
        if stats_He3[fit][0]<He3_x2_cut:
            plt.plot(r, rho_ch(r,Qich_He3[fit],R_He3[fit],2), color='red', alpha=0.2)
            Ri = R_He3[fit]
            Qich = Qich_He3[fit]
            #print('3He Charge Distribution Integral(0,10) =',quad(rho_ch_int, 0, 5, args=(Qich_He3[fit],R_He3[fit],2)))

#Dummy plot just to add a label for the ensemble fits.
plt.plot(0, 0, color='red',label='Ensemble Fits')

#Plot Amroun charge FF representative fit.
if show_amroun==1:
    plt.plot(r, rho_ch(r,Qich_He3_Amroun,R_He3_Amroun,2), color='blue',label='Amroun Fit')
    print('Amroun 3He Charge Distribution Integral(0,10) =',quad(rho_ch_int, 0, 5, args=(Qich_He3_Amroun,R_He3_Amroun,2)))

ax.legend(loc='upper right')
plt.show()


#Plot distribution of charge radii.
#User defined Gaussian fit.
fig, ax = plt.subplots(figsize=(12,6))
ax.set_ylabel('Occurrences',fontsize=16)
ax.set_xlabel('RMS Radius (fm)',fontsize=16)

xmin = 1.85 #thesis -> 1.89. sum q1 = 1 -> 1.85
xmax = 1.875 #thesis -> 1.915.sum q1 = 1 -> 1.875
nbins = 50

#Define range and number of bins.
bins = np.linspace(xmin, xmax, nbins)

#Create a histogram and fill the bins with the radii data.
data_entries_1, bins_1 = np.histogram(radii_3He, bins=bins)

#Define where the bin centers are for fitting.
binscenters = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])

#Fit the histogram with the Gaussian fit and some starting parameters.
popt, pcov = curve_fit(gaus_fit, xdata=binscenters, ydata=data_entries_1, p0=[1, 1.9, 0.1])
print('3He popt =',popt)
#print('pcov =',pcov)

#Add the average value and standard deviation of the charge radii to the title. 
ax.set_title('$^3$He Charge Radii: Average={:.3f} Standard Deviation={:.4f}'.format(popt[1], popt[2]),fontsize=20)

#Define enough points to make a smooth curve.
xspace = np.linspace(xmin, xmax, 100000)

#Plot the histogram.
plt.bar(binscenters, data_entries_1, width=bins[1] - bins[0], color='navy', label=r'Histogram entries')

#Plot the Gaussian fit to the histogram.
plt.plot(xspace, gaus_fit(xspace, *popt), color='darkorange', linewidth=2.5, label=r'Fitted function')

#Display both the histogram and the Gaussian fit.
plt.show()

#Plot the 3He magnetic FF Fits.
fig, ax = plt.subplots(figsize=(12,6))
ax.set_title('$^3$He Magnetic Form Factor',fontsize=20)
ax.set_ylabel('$F_{m}(Q^2)$',fontsize=16)
ax.set_xlabel('$Q^2$ (fm$^{-2}$)',fontsize=16)
ax.set_yscale('log')

min = 0
max = 60
plt.xticks(np.arange(min, max+1, 2.0))

#Define Q2eff range to plot.
Q2eff = np.linspace(0.00001,60,600)

#Plot ensemble of 3He fits surviving x^2 cut.
if show_ensemble==1:
    for fit in range(0,len(R_He3)):
        if stats_He3[fit][0]<He3_x2_cut:
            plt.plot(Q2eff, np.absolute(Fch(Q2eff,Qim_He3[fit],R_He3[fit])), color='red', alpha=0.2)

#Dummy plot just to add a label for the ensemble fits.
plt.plot(0, 0, color='red',label='Ensemble Fits')

#Plot Amroun magnetic FF representative fit.
if show_amroun==1:
    plt.plot(Q2eff, np.absolute(Fch(Q2eff,Qim_He3_Amroun,R_He3_Amroun)), color='blue',label='Amroun Fit')
    #This works to fill the error bands. Need to flip one error band array so the fill works correctly.
    x = np.append(He3_Fm_Amroun_Error_Band_Down_x,np.flip(He3_Fm_Amroun_Error_Band_Up_x))
    y = np.append(He3_Fm_Amroun_Error_Band_Down_y,np.flip(He3_Fm_Amroun_Error_Band_Up_y))
    plt.fill(x,y)

#Plot new 3He magnetic FF representative fit.
if show_rep==1:
    plt.plot(Q2eff, np.absolute(Fch(Q2eff,Qim_He3_thesis,R_He3_thesis)), color='black',label='New Representative Fit') #Plot 3He representative fit.

#Plot theory curves.
if show_theory==1:
    plt.plot(He3_fm_conv_x, He3_fm_conv_y, color='green',label='Conventional Approach Marcucci 2016')
    plt.plot(He3_fm_CST_x, He3_fm_CST_y, color='cyan',label='Covalent Spectator Theorem Marcucci 2016')
    plt.plot(He3_Fm_XEFT500_x, He3_Fm_XEFT500_y, color='m',label='$\chi$EFT500 Marcucci 2016')
    plt.plot(He3_Fm_XEFT600_x, He3_Fm_XEFT600_y, color='brown',label='$\chi$EFT600 Marcucci 2016')

ax.legend(loc='upper right')

plt.show()


#Plot the 3H charge FF Fits.
ngaus = 8
fig, ax = plt.subplots(figsize=(12,6))
ax.set_title('$^3$H Charge Form Factor',fontsize=20)
ax.set_ylabel('$F_{ch}(Q^2)$',fontsize=16)
ax.set_xlabel('$Q^2$ (fm$^{-2}$)',fontsize=16)
ax.set_yscale('log')

min = 0
max = 60
plt.xticks(np.arange(min, max+1, 2.0))

#Define Q2eff range to plot.
Q2eff = np.linspace(0.00001,60,600)

#Plot ensemble of 3H fits surviving x^2 cut.
if show_ensemble==1:
    for fit in range(0,len(R_H3)):
        if stats_H3[fit][0]<H3_x2_cut:
            plt.plot(Q2eff, np.absolute(Fch(Q2eff,Qich_H3[fit],R_H3[fit])), color='red', alpha=0.2)
            Ri = R_H3[fit]
            Qich = Qich_H3[fit]
            d = derivative(Fch_deriv, 0.0015, dx=1e-5)
            radius = rms_radius(d)
            radii_3H.append(radius)
            #print ('Fch_deriv =',d,'radius =',radius)

radii_3H = np.array(radii_3H)

#Dummy plot just to add a label for the ensemble fits.
plt.plot(0, 0, color='red',label='Ensemble Fits')

#Plot Amroun charge FF representative fit.
if show_amroun==1:
    plt.plot(Q2eff, np.absolute(Fch(Q2eff,Qich_H3_Amroun,R_H3_Amroun)), color='blue',label='Amroun Fit')
    #This works to fill the error bands. Need to flip one error band array so the fill works correctly.
    x = np.append(H3_Fch_Amroun_Error_Band_Down_x,np.flip(H3_Fch_Amroun_Error_Band_Up_x))
    y = np.append(H3_Fch_Amroun_Error_Band_Down_y,np.flip(H3_Fch_Amroun_Error_Band_Up_y))
    plt.fill(x,y)

#Plot new 3H charge FF representative fit.
if show_rep==1:
    plt.plot(Q2eff, np.absolute(Fch(Q2eff,Qich_H3_thesis,R_H3_thesis)), color='black',label='New Representative Fit') #Plot 3H representative fit.

#Plot theory curves.
if show_theory==1:
    plt.plot(H3_fch_conv_x, H3_fch_conv_y, color='green',label='Conventional Approach Marcucci 2016')
    plt.plot(H3_fch_CST_x, H3_fch_CST_y, color='cyan',label='Covalent Spectator Theorem Marcucci 2016')
    plt.plot(H3_Fch_XEFT500_x, H3_Fch_XEFT500_y, color='m',label='$\chi$EFT500 Marcucci 2016')
    plt.plot(H3_Fch_XEFT600_x, H3_Fch_XEFT600_y, color='brown',label='$\chi$EFT600 Marcucci 2016')

ax.legend(loc='upper right')

plt.show()

#Plot the 3H charge distributions.
fig, ax = plt.subplots(figsize=(12,6))
ax.set_title('$^3$H Charge Density',fontsize=20)
ax.set_ylabel(r'$\rho(r)$ (e/fm$^3$)',fontsize=16)
ax.set_xlabel('Radius (fm)',fontsize=16)
#ax.set_yscale('log')

#Define radii range to plot.
r = np.linspace(0.00001,5,5000)

if show_ensemble==1:
    for fit in range(0,len(R_H3)):
        if stats_H3[fit][0]<H3_x2_cut:
            plt.plot(r, rho_ch(r,Qich_H3[fit],R_H3[fit],1), color='red', alpha=0.2)
            Ri = R_H3[fit]
            Qich = Qich_H3[fit]
            #print('3H Charge Distribution Integral(0,10) =',quad(rho_ch_int, 0, 5, args=(Qich_H3[fit],R_H3[fit],1)))

#Dummy plot just to add a label for the ensemble fits.
plt.plot(0, 0, color='red',label='Ensemble Fits')

#Plot Amroun charge FF representative fit.
if show_amroun==1:
    plt.plot(r, rho_ch(r,Qich_H3_Amroun,R_H3_Amroun,1), color='blue',label='Amroun Fit')
    print('Amroun 3H Charge Distribution Integral(0,10) =',quad(rho_ch_int, 0, 5, args=(Qich_H3_Amroun,R_H3_Amroun,1)))

ax.legend(loc='upper right')
plt.show()

#Plot distribution of charge radii.
#User defined Gaussian fit.
fig, ax = plt.subplots(figsize=(12,6))
ax.set_ylabel('Occurrences',fontsize=16)
ax.set_xlabel('RMS Radius (fm)',fontsize=16)

xmin = 1.65 #thesis -> 1.97. sum q1 = 1 -> 1.65
xmax = 1.75 #thesis -> 2.06. sum q1 = 1 -> 1.75
nbins = 50

#Define range and number of bins.
bins = np.linspace(xmin, xmax, nbins)

#Create a histogram and fill the bins with the radii data.
data_entries_1, bins_1 = np.histogram(radii_3H, bins=bins)

#Define where the bin centers are for fitting.
binscenters = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])

#Fit the histogram with the Gaussian fit and some starting parameters.
popt, pcov = curve_fit(gaus_fit, xdata=binscenters, ydata=data_entries_1, p0=[1, 1.9, 0.1])
print('3H popt =',popt)
#print('pcov =',pcov)

#Add the average value and standard deviation of the charge radii to the title. 
ax.set_title('$^3$H Charge Radii: Average={:.3f} Standard Deviation={:.4f}'.format(popt[1], popt[2]),fontsize=20)

#Define enough points to make a smooth curve.
xspace = np.linspace(xmin, xmax, 100000)

#Plot the histogram.
plt.bar(binscenters, data_entries_1, width=bins[1] - bins[0], color='navy', label=r'Histogram entries')

#Plot the Gaussian fit to the histogram.
plt.plot(xspace, gaus_fit(xspace, *popt), color='darkorange', linewidth=2.5, label=r'Fitted function')

#Display both the histogram and the Gaussian fit.
plt.show()

#Plot the 3H magnetic FF Fits.
fig, ax = plt.subplots(figsize=(12,6))
ax.set_title('$^3$H Magnetic Form Factor',fontsize=20)
ax.set_ylabel('$F_{m}(Q^2)$',fontsize=16)
ax.set_xlabel('$Q^2$ (fm$^{-2}$)',fontsize=16)
ax.set_yscale('log')

min = 0
max = 60
plt.xticks(np.arange(min, max+1, 2.0))

#Define Q2eff range to plot.
Q2eff = np.linspace(0.00001,60,600)

#Plot ensemble of 3He fits surviving x^2 cut.
if show_ensemble==1:
    for fit in range(0,len(R_H3)):
        if stats_H3[fit][0]<H3_x2_cut:
            plt.plot(Q2eff, np.absolute(Fch(Q2eff,Qim_H3[fit],R_H3[fit])), color='red', alpha=0.2)

#Dummy plot just to add a label for the ensemble fits.
plt.plot(0, 0, color='red',label='Ensemble Fits')

#Plot Amroun magnetic FF representative fit.
if show_amroun==1:
    plt.plot(Q2eff, np.absolute(Fch(Q2eff,Qim_H3_Amroun,R_H3_Amroun)), color='blue',label='Amroun Fit')
    #This works to fill the error bands. Need to flip one error band array so the fill works correctly.
    x = np.append(H3_Fm_Amroun_Error_Band_Down_x,np.flip(H3_Fm_Amroun_Error_Band_Up_x))
    y = np.append(H3_Fm_Amroun_Error_Band_Down_y,np.flip(H3_Fm_Amroun_Error_Band_Up_y))
    plt.fill(x,y)

#Plot new 3H magnetic FF representative fit.
if show_rep==1:
    plt.plot(Q2eff, np.absolute(Fch(Q2eff,Qim_H3_thesis,R_H3_thesis)), color='black',label='New Representative Fit') #Plot 3H representative fit.

#Plot theory curves.
if show_theory==1:
    plt.plot(H3_fm_conv_x, H3_fm_conv_y, color='green',label='Conventional Approach Marcucci 2016')
    plt.plot(H3_fm_CST_x, H3_fm_CST_y, color='cyan',label='Covalent Spectator Theorem Marcucci 2016')
    plt.plot(H3_Fm_XEFT500_x, H3_Fm_XEFT500_y, color='m',label='$\chi$EFT500 Marcucci 2016')
    plt.plot(H3_Fm_XEFT600_x, H3_Fm_XEFT600_y, color='brown',label='$\chi$EFT600 Marcucci 2016')

ax.legend(loc='upper right')

plt.show()


