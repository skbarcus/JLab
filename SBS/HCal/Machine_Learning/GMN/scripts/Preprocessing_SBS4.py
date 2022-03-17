import numpy as np
import tensorflow as tf
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split

from matplotlib.ticker import MaxNLocator #Able to force plots to use integer ticks only. 
from keras.utils import np_utils #Utilities to transform data.
from sklearn.preprocessing import normalize

import time
start = time.time() #Start a timer.

save_data = 1
data = '/home/skbarcus/JLab/SBS/HCal/Machine_Learning/GMN/data/'

#Read in the data line by line.
with open('/home/skbarcus/JLab/SBS/HCal/Machine_Learning/GMN/data/elastic_training_data_SBS4.txt') as f:
    lines = f.readlines()

#Remove first line with column labels.
del lines[0]

#Check data shape.
print(lines[0])
print(lines[1])
print(lines[2])

#Create array to store data.
SBS4_Data = []

#Read each line and split by the spaces and store as array.
for line in lines:
    event = line.split()
    SBS4_Data.append(event)
    #print(event)

#Turn data into numpy array.
SBS4_Data = np.array(SBS4_Data)

#Examine data shape and check output.
print('SBS4_Data.shape = ',SBS4_Data.shape)
print('SBS4_Data = ',SBS4_Data)

#Split array into targets (elastic or not elastic) and training data.
np.hsplit(SBS4_Data,np.array([1,6]))

#Examine data shape and check output.
print('SBS4_Data.shape = ',SBS4_Data.shape)
print('SBS4_Data = ',SBS4_Data)

#Split off first column into its own array for targets and last 6 cols as training data.
SBS4_Data = np.hsplit(SBS4_Data,np.array([1]))

target_data = np.array(SBS4_Data[0].astype('float'))
training_data = np.array(SBS4_Data[1].astype('float'))

#print(len(SBS4_Data))
#for i in range(0,len(SBS4_Data)):
#    target_data[i] = SBS4_Data[i][0]
#    training_data[i] = 

#print('training_data[0][0] = ',training_data[0][0])
#print('training_data[1][0] = ',training_data[1][0])
#print('training_data[2][0] = ',training_data[2][0])
#print('training_data[9][0] = ',training_data[9][0])

#Remove events where HCal cluster energy is zero. Bad idea. Very very slow. Clean data earlier.
#for i in range(0,len(training_data)):
#    if i%50000==0:
#        print(i)
#    if training_data[i][0].astype('float') == 0:
#        #print('Hello')
#        target_data = np.delete(target_data,i,0)
#        training_data = np.delete(training_data,i,0)

training_data_temp = []
target_data_temp = []

#Find min and max values for each variable in the training data. Start by setting vectors equal to first event.
max_training = []
min_training = []
for i in range(0,6):
    max_training.append(training_data[0][i])
    min_training.append(training_data[0][i])

print('Initial max_training =',max_training)
print('Initial min_training =',min_training)

#Remove events where HCal cluster energy is zero or tdc time is -1000.
for i in range(0,len(training_data)):
    if i%50000==0:
        print(i)
    if training_data[i][0].astype('float') != 0 and training_data[i][1].astype('float') != -1000:
        training_data_temp.append(training_data[i])
        target_data_temp.append(target_data[i])
#        max_training_temp = []
#        for j in range(0,len(max_training)):
#            if training_data[i][j]>max_training[j]:
#                max_training[j]=training_data[i][j]
#            else:
#                max_training[j]=max_training[j]
#            if training_data[i][j]<min_training[j]:
#                min_training[j]=training_data[i][j]
#            else:
#                min_training[j]=min_training[j]
        
training_data_temp = np.array(training_data_temp)
target_data_temp = np.array(target_data_temp)

print('training_data_temp.shape = ',training_data_temp.shape)
print('training_data_temp = ',training_data_temp)

print('target_data_temp.shape = ',target_data_temp.shape)
print('target_data_temp = ',target_data_temp)

#Replace old arrays with the cleaned up temp arrays.
target_data = target_data_temp
training_data = training_data_temp

numpy_max = np.amax(training_data,axis=0)
print('Numpy Array Maxes',numpy_max)
numpy_min = np.amin(training_data,axis=0)
print('Numpy Array Mins',numpy_min)

print('target_data.shape = ',target_data.shape)
print('target_data = ',target_data)
#target_data = normalize(target_data)
#print('Normalized target_data = ',target_data)

print('training_data.shape = ',training_data.shape)
print('training_data = ',training_data)

#Normalize the data so it's useful for neural network learning.
#Moved this fucntion to the notebook soI have raw data to plot and analyze.
#training_data = normalize(training_data,axis=0)
#print('Normalized training_data = ',training_data)


#print('max_training = ',max_training)
#print('min_training = ',min_training)

if save_data == 1:
    np.save(data + "target_data_SBS4", target_data)
    np.save(data + "training_data_SBS4", training_data)

print("The script took %.2f seconds (%.2f minutes or %.2f hours) to run." % (time.time() - start, (time.time() - start)/60.,(time.time() - start)/60./60.)) #Print time to run.
