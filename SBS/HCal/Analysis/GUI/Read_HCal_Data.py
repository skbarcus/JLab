from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
from mpl_toolkits.axes_grid1.mpl_axes import Axes
import numpy as np
import ROOT
import sys

from keras.utils import np_utils #Utilities to transform data.

from sklearn.preprocessing import normalize

save_data = 1
read_all = 0
loop = 100000
test_evt = 2#31823
nsamps = 20

data = "/home/skbarcus/Machine_Learning/Tutorials/Data/"

infile = sys.argv[1]

print("Reading from", infile)

inFile = ROOT.TFile.Open(infile," READ ")

tree = inFile.Get("T")

h1 = ROOT.TH1D("h1","h1 test",20,0,50000)

print(tree.GetEntries())
if read_all == 1:
    loop = tree.GetEntries()

adc_samps_vals = []
adc_vals = []
tdc_vals = []
hit_vals = []

#Loop over entries and fill histo.
for entryNum in range (0, loop):
    tree.GetEntry(entryNum)
    adc_samps_evt = getattr(tree,"sbs.hcal.samps")
    adc_samps_vals.append([])
    #print('adc_samps_evt = ', adc_samps_evt)
    #print('adc_samps_evt[0] = ', adc_samps_evt[0])
    adc_evt = getattr(tree,"sbs.hcal.a")
    #print('adc_evt = ', adc_evt)
    adc_vals.append([])
    tdc_evt = getattr(tree,"sbs.hcal.tdc")
    tdc_vals.append([])
    #print(type(adc_evt))
    #print(adc_evt[0]," ",adc_evt[1])
    hit_flag = 0
    #hit.append([])
    #Loop over PMTs.
    for samp in range(0,len(adc_samps_evt)):
        adc_samps_vals[entryNum].append(adc_samps_evt[samp])
    for pmt in range(0,len(adc_evt)):
        adc_vals[entryNum].append(adc_evt[pmt])
        tdc_vals[entryNum].append(tdc_evt[pmt])
        h1.Fill(adc_evt[pmt])  
        if adc_evt[pmt] > 7000 and adc_evt[pmt] < 25000:
            hit_flag = 1
    if hit_flag == 1:
        hit_vals.append(1)
    else:
        hit_vals.append(0)
            
#print(adc_vals[0:3])
#print(hit_vals[0:5])

#Convert lists into numpy arrays.
adc_samps_arr = np.array(adc_samps_vals)
print('adc_samps_arr shape = ',adc_samps_arr.shape)
#print(adc_samps_arr[0][:40])

adc_arr = np.array(adc_vals)
print('adc_arr shape = ',adc_arr.shape)

tdc_arr = np.array(tdc_vals)
print('tdc_arr shape = ',tdc_arr.shape)

hit_arr = np.array(hit_vals)
print('hit_arr shape = ',hit_arr.shape)

#hit_arr = np_utils.to_categorical(hit_arr, 2)
#print('hit_arr shape = ',hit_arr.shape)

#print('Unormalized adc_arr = ',adc_arr)
#print(adc_arr[0][0])
#print(adc_arr[0][1])
#print(adc_arr[0,0])

# Plot the two pixel arrays side by side.
"""
fig, ax1 = plt.subplots()
#ax.plot(adc_arr[2])#Plots the fadc integral for each of the 144 PMTs for the adc_arr[event]. x-axis = PMT #. y-axis = fadc integral.
#ax.plot(adc_arr)#Plots 144 plots (one per PMT) of their fadc integral value for each event. x-axis = event #. y-axis fadc integral.

ax1.plot(adc_arr[:,0:3])#Plots fadc integral for all events for PMTs 0, 1, 2.
plt.show()
"""

#Set the fadc integrals of the two reference channels to 5000 for ease of plotting.
adc_arr[:,142:144] = 5000

#test_pmt = 68
# print('unormalized test_pmt = ',adc_arr[test_evt][test_pmt])
# norm_test = normalize(adc_arr)
# print('normalized test_pmt = ',norm_test[test_evt][test_pmt])
# print(norm_test)

#adc_arr = normalize(adc_arr)
#print('Normalized adc_arr = ', adc_arr)

#Plot all fADC integrals as a histogram. 
"""
Do NOT use this. If plot more than a few 10000s events will eat up all of laptop's RAM.
cfadc_int = ROOT.TCanvas("cfadc_int","fADC Integrals for All Events",900,700)
hfadc_int = ROOT.TH1F("hfadc_int","fADC Integrals for All Events",1000,0,165000)
for i in range (0,loop):
    for j in range (0,len(adc_arr[1])):
        hfadc_int.Fill(adc_arr[i][j])
cfadc_int.cd(1)
hfadc_int.Draw("hist")
cfadc_int.Update()
input("Please press enter to close images.")  
"""

#Plot fADC integrals for a single event. 
fig, ax2 = plt.subplots()
adc_int_event = ax2.imshow(adc_arr[test_evt].reshape(12,12))
fig.colorbar(adc_int_event, ax=ax2)
#ax2.plot(tdc_arr[:,0])
plt.show()

#Plot all fADC samples for a single event. 
fig, ax3 = plt.subplots(figsize=(10,10))
#fig, ax3 = plt.subplots()
adc_samps_event = ax3.imshow(adc_samps_arr[test_evt].reshape(12,240),aspect='auto')
fig.colorbar(adc_samps_event, ax=ax3)
#ax2.plot(tdc_arr[:,0])
plt.show()

xaxis = []
for i in range(nsamps):
    xaxis.append(i)
#print('xaxis = ',xaxis)
yaxis = []
for i in range(nsamps):
    yaxis.append(adc_samps_arr[test_evt][i])
#print('yaxis = ',yaxis)

#Plot 1D ROOT histogram of a single PMT fADC channel.
h1dfadc = ROOT.TH1F("hfadc","hfadc test",20,0,20)
for i in range (0,nsamps):
    h1dfadc.Fill(xaxis[i],yaxis[i])

h1dfadc.Draw("hist")
input("Please press enter to close images.")  

#Plot 2D ROOT histogram of a single PMT fADC channel. Less useful.
"""
h2dfadc = ROOT.TH2F("h2dfadc","h2dfadc test",20,0,20,1000,0,8200)
for samps in range (0,nsamps):
    h2dfadc.Fill(samps,adc_samps_arr[test_evt][samps])
    print(samps,"   ",adc_samps_arr[test_evt][samps])

h2dfadc.Draw("box")
input("Please press enter to close images.")  
"""

#Plot 1D Matplot histogram of a single PMT fADC channel.
"""
fig, ax = plt.subplots(figsize=(10,10))
hist1 = ax.hist(xaxis,bins=nsamps,weights=yaxis)
plt.show()
"""

#Plot 2D Matplot histogram of a single PMT fADC channel. Less useful.
"""
fig, ax = plt.subplots(figsize=(10,10))
hist = ax.hist2d(xaxis, yaxis)
plt.show()
"""

#Plot all fADC sample histos for a single event.
#Create list to store a root hist for each PMT. 
histos = []
#Create root hist for each PMT's fadc samples for a single event.
for i in range(0,adc_arr.shape[1]):
    histos.append(ROOT.TH1F('hfadc%d' % (i),'fADC Samples for PMT %d' % (i),20,0,20))
#print('histos = ',histos)
#Fill the root histograms with the fADC samples for each PMT.
for i in range(0,len(histos)):
    yaxis = []
    for j in range(nsamps):
        yaxis.append(adc_samps_arr[test_evt][i*nsamps+j])
    for j in range(nsamps):
        histos[i].Fill(xaxis[j],yaxis[j])

ctest = ROOT.TCanvas("ctest","ctest",900,700)
#ROOT.gStyle.SetOptStat(0)#Might be an issue.
ctest.Divide(12,12,0,0)#(cols,rows,0,0)#Fails above (3,3,0,0). (2,3) and (3,2) work.
"""
ctest.cd(1)
#ctest.Update()
histos[0].Draw("hist")
ctest.Update()#This is important for displaying the histos.
"""

for i in range(0,len(histos)):
    ctest.cd(i+1)
    histos[i].Draw("hist")
ctest.Update()
#histos[0].Draw("hist") 
#fig, ax = plt.subplots(12,12,figsize=(10,10))
#for i in range(0,len(histos)):
    
input("Please press enter to close images.")
"""
#Fills a single histogram.
print('histos[0] = ',histos[0])
for i in range(0,nsamps):
    histos[0].Fill(xaxis[i],yaxis[i])
histos[0].Draw("hist") 
input("Please press enter to close images.")
"""
#Normalize the arrays so they're compatible with NNs.
adc_arr = normalize(adc_arr)
#print('Normalized adc_arr = ', adc_arr)
adc_samps_arr = normalize(adc_samps_arr)
#print('Normalized adc_samps_arr = ', adc_samps_arr)

if save_data == 1:
    np.save(data + "adc_integrals_run820", adc_arr)
    np.save(data + "adc_samps_run820", adc_samps_arr)
    np.save(data + "tdc_run820", tdc_arr)
    np.save(data + "hits_run820", hit_arr)

