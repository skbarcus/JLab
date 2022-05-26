#Script to create a bootstrapped dataset from an original. Samples the original randomly with replacement to create a new dataset with the same number of data points. 

import numpy as np
import matplotlib.pyplot as plt
#Generate random integer values
from random import seed
from random import random
from random import randint

#Read in data.
with open('/home/skbarcus/JLab/SOG/3He_Data.txt') as f:
#with open('/home/skbarcus/JLab/SOG/3H_Data_Thesis.txt') as f:
    lines = f.readlines()

#Remove lines with column labels.
del lines[0]
del lines[0]

#Create arrays.
Raw_Data = []          #Array to store raw data.
Bootstrap_Data = []    #Array to store bootstrapped data.

#Read each line and split by the spaces and store as array.
#['Energy (GeV)', 'Theta (Degrees)', 'Sigma Experimental', 'Uncertainties', 'Dataset']
for line in lines:
    event = line.split()
    Raw_Data.append(event)

#Turn data into numpy array and swap char to float.
Raw_Data = np.array(Raw_Data)
Raw_Data = np.array(Raw_Data.astype('float'))

#Examine data shape and check output.
print('Raw_Data.shape = ',Raw_Data.shape)
print('Raw_Data[0] = ',Raw_Data[0])

#Create text file for bootstrapped data.
with open('He3_Bootstrapped_Data8888.txt', 'w') as f:
    #Add some header lines.
    f.write('Bootstrapped Data\n')
    f.write('Energy (GeV)   Theta (Degrees)   Sigma Experimental   Uncertainties   Dataset\n')

    #Sample randomly from raw data with replacement and build a new array of bootstrapped data points.
    for i in range(0,len(Raw_Data)):
        rand_idx = randint(0, len(Raw_Data)-1)
        #print('New point',i,': uses original index',rand_idx)
        Bootstrap_Data.append(Raw_Data[rand_idx])
        #Write bootstrapped data to text file.
        for j in range(0,4):
            f.write(Bootstrap_Data[i][j].astype('str'))
            f.write('   ')
        f.write('\n')

#Turn data into numpy array and swap char to float.
Bootstrap_Data = np.array(Bootstrap_Data)

#Examine data shape and check output.
print('Bootstrap_Data.shape = ',Bootstrap_Data.shape)
print('Bootstrap_Data[0] = ',Bootstrap_Data[0])

