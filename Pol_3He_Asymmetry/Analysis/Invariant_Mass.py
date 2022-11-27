from matplotlib import pyplot as plt
import numpy as np
import math
import time
from matplotlib.transforms import Transform
from matplotlib.ticker import(AutoLocator, AutoMinorLocator)
from collections import OrderedDict

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
MtHe3 = 3.0160293*0.9315        #3He mass GeV.
MtH = 0.938                     #H mass GeV.

def W(theta,E0,Ef):
    theta = theta*deg2rad
    #E0 = E
    #Ef = E0/(1.0+2.0*E0*pow(np.sin(theta/2.0),2.0)/MtHe3)
    Q2 = 4*E0*Ef*pow(np.sin(theta/2),2)
    W = pow( abs(pow(MtHe3,2) + 2*MtHe3*(E0-Ef)-Q2) ,0.5)
    print("Ef =",Ef)
    print("Q^2 =",Q2,"GeV^2 = ",Q2*25.7,"fm^(-2).")
    print("W =",W)
    return W

W(13,2.18,2)
