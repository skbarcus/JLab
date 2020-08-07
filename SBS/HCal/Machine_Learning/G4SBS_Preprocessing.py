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
read_all = 1
loop = 20000
data = "/home/skbarcus/JLab/SBS/HCal/Analysis/Simulation/G4SBS/rootfiles/"
npmts = 288       #Number of PMTs/Scintillators in G4SBS HCal simulation.

def plot_learning_curve(history):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    ax[0].plot(history['loss'], label='training')
    ax[0].plot(history['val_loss'], label='validation')
    ax[0].set_title("Model losses")
    ax[0].set_xlabel("Epoch")
    ax[0].set_ylabel("Loss")
    ax[0].set_yscale('log')
    ax[0].xaxis.set_major_locator(MaxNLocator(integer=True))
    ax[0].legend()
    
    ax[1].plot(history['accuracy'], label='training')
    ax[1].plot(history['val_accuracy'], label='validation')
    ax[1].set_title("Model accuracy")
    ax[1].set_xlabel("Epoch")
    ax[1].set_ylabel("Accuracy")
    ax[1].xaxis.set_major_locator(MaxNLocator(integer=True))
    ax[1].legend()
    plt.show()

data = "/home/skbarcus/JLab/SBS/HCal/Analysis/Simulation/G4SBS/rootfiles/"

fnucl = np.load(data+"g4sbs_kin7_fnucl.npy")
front_nhits = np.load(data+"g4sbs_kin7_front_nhits.npy")
nhits = np.load(data+"g4sbs_kin7_nhits.npy")
pmt = np.load(data+"g4sbs_kin7_pmt.npy", allow_pickle=True)
pmt_edep_data = np.load(data+"g4sbs_kin7_pmt_edep.npy", allow_pickle=True)
tavg_data = np.load(data+"g4sbs_kin7_tavg.npy", allow_pickle=True)

print("fnucl shape = ",fnucl.shape)
print("front_nhits shape = ",front_nhits.shape)
print("nhits shape = ",nhits.shape)
print("pmt shape = ",pmt.shape)
print("pmt_edep_data shape = ",pmt_edep_data.shape)
print("tavg_data shape = ",tavg_data.shape)
#print(fnucl)
#print(tavg)
#print(pmt)
#print(type(pmt[0]))

pmt_edep = []
tavg = []

if read_all == 1:
    loop = len(fnucl)

for i in range(0,loop):
    pmt_edep.append([])
    tavg.append([])
    idx = 0     #Track the index of the stored arrays that aren't well shaped.
    for j in range(0,npmts):
        if j in pmt[i]:
            #print(j)
            pmt_edep[i].append(pmt_edep_data[i][idx])
            tavg[i].append(tavg_data[i][idx])
            idx = idx +1
        else:
            pmt_edep[i].append(0)
            tavg[i].append(0)
    if i%10000==0:
        print("Read in %d events." % i)

#Normalize the arrays so they're compatible with NNs.
pmt_edep = np.array(pmt_edep)
#print(pmt_edep)
pmt_edep = normalize(pmt_edep)
print("pmt_edep shape = ", pmt_edep.shape)
#print(pmt_edep_data[0])
#print(pmt_edep[0])

tavg = np.array(tavg)
#print(tavg)
tavg = normalize(tavg)
print("tavg shape = ", tavg.shape)
#print(tavg_data[0])
#print(tavg[0])

#Fill all arrays with the scintillators/pmts that saw no signal so the array shapes are (***,288).
"""
for i in range(0,len(fnucl)):
    for j in range(0,len(pmt[i])):
        for k in range(0,npmts):
            print(pmt[i][k])
            #if j==pmt[i][j]:
"""
"""
a = np.array([[1,2,3],[4,5,6]]) 

print('First array:') 
print(a) 
print('\n')  

print('Append elements to array:') 
print(np.append(a, [7,8,9])) 
print('\n')  

print('Append elements along axis 0:') 
print(np.append(a, [[7,8,9]],axis = 0)) 
print('\n')  

print('Append elements along axis 1:') 
print(np.append(a, [[5,5,5],[7,8,9]],axis = 1))
"""
if save_data == 1:
     np.save(data + "g4sbs_kin7_pmt_edep_preprocessed", pmt_edep)
     np.save(data + "g4sbs_kin7_tavg_preprocessed", tavg)

print("The script took %.2f seconds (%.2f minutes or %.2f hours) to run." % (time.time() - start, (time.time() - start)/60.,(time.time() - start)/60./60.)) #Print time to run.
