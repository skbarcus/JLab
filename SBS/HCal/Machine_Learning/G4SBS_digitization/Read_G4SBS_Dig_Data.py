from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
from mpl_toolkits.axes_grid1.mpl_axes import Axes
# generate random Gaussian values
from random import seed
from random import gauss
import numpy as np
import ROOT
#from ROOT import TChain, TSelector, TTree
import sys

#from keras.utils import np_utils #Utilities to transform data.
#from sklearn.preprocessing import normalize
#from root_numpy import root2array, testdata
#root2array(testdata.get_filepath('single1.root'))[:20]

import time
start = time.time() #Start a timer.

save_data = 0
read_all = 1
loop = 30000
test_evt = 2#31823
nsamps = 20
npmts = 288

data = "/home/skbarcus/JLab/SBS/HCal/Machine_Learning/rootfiles/"
#infile = sys.argv[1] #Use file specified when code is run. 

"""
chain = ROOT.TChain("T")
for i in range(0,20):
    rootfile = "simdigtest_"+("%d" % i)+".root"
    chain.AddFile(data+rootfile)
    print("Adding rootfile ",data+rootfile," to chain.")
"""

chain = ROOT.TChain("T")
rootfile = "simdigtest_0.root"
chain.AddFile(data+rootfile)
print("Opening rootfile ",data+rootfile,".")

print(chain.GetEntries())
if read_all == 1:
    loop = chain.GetEntries()

print(rootfile,"has ",loop,"entries.")

fnucl = []         #0 = neutron, 1 = proton.
front_nhits = []   #Hits on the front panel of HCal.
nhits = []         #Scintillator hits.
pmt = []           #PMT module with a hit.
pmt_edep = []      #Energy deposited in a single PMT (scintillator).
tavg = []          #Average time scint hit occurred.
dig_adc = []       #Holds digitized fADC values in an array of size dig_adc[nsamps] with dig_adc[0] holding all ADC values for the first fADC bin in an event.
pmt_adc = []       #Holds the fADC values for each PMT in each event. If G4SBS said PMT didn't fire, it fills the PMT with random pedestal noise.
pmt_hits = []      #1 if the PMT fired according to G4SBS and 0 if it didn't.
adc_int = []       #Array to hold all 288 fADC values for each event.
tdc = []
for i in range(0,nsamps):
    dig_adc.append([])
srau = []          #Integrals of PMT fADCs.
chan = []          #Stores the PMT channel of the hit. Values 0-287.

#Loop over entries.
loop = 1
for entryNum in range (0, loop):
    chain.GetEntry(entryNum)

    fnucl_evt = getattr(chain,"fnucl")#0 = neutron, 1 = proton.
    fnucl.append(fnucl_evt)

    chan_evt = getattr(chain,"Harm.HCal.dighit.chan")
    #print(chan_evt)
    chan_out = np.array(chan_evt)
    chan.append(chan_out)
    #print(chan_out)

    tdc_evt = getattr(chain,"Harm.HCal.dighit.tdc")
    #print(len(tdc_evt),tdc_evt)
    tdc_out = np.array(tdc_evt)
    tdc.append(tdc_out)
    #print(tdc_out)

    for samp in range(0,nsamps):
        current_branch = "Harm.HCal.dighit.adc_"+str(samp)
        adc_evt = getattr(chain,current_branch)
        #print("adc_evt",adc_evt)
        adc_evt = np.array(adc_evt)
        dig_adc[samp].append(adc_evt)
        #current_samp = "dig_adc_"+str(samp)
        #current_samp.append(adc)
        #dig_adc_0_evt = getattr(chain,"Harm.HCal.dighit.adc_0")
        #adc = np.array(dig_adc_0_evt)
        #dig_adc_0.append(adc)

    #adc = np.array(dig_adc_0_evt)
    #print(adc)
    if entryNum%5000==0:
        print('Processing digitized rootfile is '+str((entryNum/loop)*100.)+'% complete.')
    elif entryNum==(loop-1):
        print('Processing digitized rootfile is 100% complete.')
        print("This took %.2f seconds (%.2f minutes or %.2f hours)." % (time.time() - start, (time.time() - start)/60.,(time.time() - start)/60./60.))
chan = np.array(chan)
#print(chan.shape)
#print('chan =',chan)

#Fill arrays of size pmt_hits[events][288] with 0 or 1 based on if the PMT fired.
for event in range(0,loop):
    pmt_hits_evt = []
    for pmt in range(0,npmts):
        if pmt in chan[event]:
            pmt_hits_evt.append(1)
        else:
            pmt_hits_evt.append(0)
    pmt_hits.append(pmt_hits_evt)
    if event%5000==0:
        print('Building PMT hits array is '+str((event/loop)*100.)+'% complete.')
    elif event==(loop-1):
        print('Building PMT hits array is 100% complete.')
        print("This took %.2f seconds (%.2f minutes or %.2f hours)." % (time.time() - start, (time.time() - start)/60.,(time.time() - start)/60./60.))

pmt_hits = np.array(pmt_hits)
#print('pmt_hits',pmt_hits)

#dig_adc[fADC samp][event][hit in event]
dig_adc = np.array(dig_adc)
for i in range(0,nsamps):
    dig_adc[i] = np.array(dig_adc[i])
#print("ADC 0 Values:")
#print(dig_adc.shape)
#print(dig_adc[0])

#Fill the array with fADC integrals.
#print("dig_adc.shape",dig_adc.shape)
for event in range(0,dig_adc.shape[1]):#Loop over events.
    integrals = np.zeros(dig_adc[samp][event].shape)
    #integrals.reshape()
    for samp in range(0,dig_adc.shape[0]):#Loop over all fADC samples in window.
        integrals = np.add(integrals,dig_adc[samp][event])
    srau.append(integrals)
    if event%5000==0:
        print('Building G4SBS PMT fADC integral values array is '+str((event/loop)*100.)+'% complete.')
    elif event==(loop-1):
        print('Building G4SBS PMT fADC integral values array is 100% complete.')
        print("This took %.2f seconds (%.2f minutes or %.2f hours)." % (time.time() - start, (time.time() - start)/60.,(time.time() - start)/60./60.))

srau = np.array(srau)
#print('srau =',srau)

#Fill arrays of size adc_int[events][288] with the fADC integrals for each PMT.
for event in range(0,loop):
    adc_int_evt = []
    hit_idx = 0
    for pmt in range(0,npmts):
        if pmt in chan[event]:
            #print('srau['+str(event)+']['+str(hit_idx)+'] =',srau[event][hit_idx])
            adc_int_evt.append(srau[event][hit_idx])
            hit_idx = hit_idx+1
        else:
            adc_int_evt.append(round(gauss(300, 1.6)*nsamps))
    adc_int.append(adc_int_evt)
    if event%5000==0:
        print('Building all PMT fADC integral values array is '+str((event/loop)*100.)+'% complete.')
    elif event==(loop-1):
        print('Building all PMT fADC integral values array is 100% complete.')
        print("This took %.2f seconds (%.2f minutes or %.2f hours)." % (time.time() - start, (time.time() - start)/60.,(time.time() - start)/60./60.))

adc_int = np.array(adc_int)
#print('adc_int',adc_int)

#for i in range(0,nsamps):
    #print(dig_adc[i][0][3])

#print("Full digitized ADC:",dig_adc)
#print(len(fnucl))
#print(fnucl)

def plot_pmt(event,hit):
    #dig_adc[fADC samp][event][hit in event]
    pmt = chan[event][hit]
    print("PMT",pmt)
    xaxis = []
    for i in range(0,nsamps):
        xaxis.append(i)
    print("xaxis =",xaxis)

    yaxis = []
    for i in range(0,nsamps):
        yaxis.append(dig_adc[i][event][hit])
    print("yaxis =",yaxis)

    #Plot 1D ROOT histogram of a single PMT fADC channel.
    hfadc = ROOT.TH1F("hfadc","hfadc",20,0,20)
    for i in range (0,nsamps):
        hfadc.Fill(xaxis[i],yaxis[i])

    hfadc.Draw("hist")
    input("Please press enter to close images.")  

#Plot all fADC sample histos for a single event.
def plot_event(event):
    #Create list to store a root hist for each PMT.
    histos = []
    ped_histos = []
    #Create and fill root histos for each PMT's fadc samples for a single hit event. Maybe add cuts on eng later.
    #print("len(dig_adc[0][event]) =",len(dig_adc[0][event]))
    pmts = []
    hit_idx = 0 
    #for hit in range(0,len(dig_adc[0][event])):
    for pmt in range(0,npmts):
        #pmt = chan[event][hit]
        histos.append(ROOT.TH1F('hfadc%d' % (pmt+1),'PMT %d' % (pmt+1),20,0,20))
        #histos[pmt].SetTitleSize(10.,"t");
        ROOT.gStyle.SetTitleSize(1,"t");
        ROOT.gStyle.SetTitleY(1.5)
        histos[pmt].SetStats(0);
        histos[pmt].GetYaxis().SetLabelSize(0.15);
        histos[pmt].GetYaxis().SetLabelOffset(-0.2);
        histos[pmt].GetYaxis().SetNdivisions(5);
        histos[pmt].SetLineWidth(2)
        if pmt in chan[event]:
            #print(pmt+1)
            xaxis = []
            yaxis = []
            #Color the PMTs with hits red.
            histos[pmt].SetLineColor(2)
            for samp in range(0,nsamps):
                xaxis.append(samp)
                yaxis.append(dig_adc[samp][event][hit_idx])
                #print("xaxis =",xaxis)
                #print("yaxis =",yaxis)

            #Fill the individual fADC hit histogram.
            for samp in range(0,nsamps):
                histos[pmt].Fill(xaxis[samp],yaxis[samp])
            hit_idx = hit_idx+1

        else:#If no hit, fill the histo with pedestal.
            xaxis = []
            yaxis = []
            for samp in range(0,nsamps):
                xaxis.append(samp)
                pedestal = round(gauss(300, 1.6))
                yaxis.append(pedestal)
            for samp in range(0,nsamps):
                histos[pmt].Fill(xaxis[samp],yaxis[samp])

    c1 = ROOT.TCanvas("c1","Subassembly 1",1200,700)
    c1.Divide(12,6)
    for pmt in range(0,72):
        c1.cd(pmt+1)
        histos[pmt].Draw("hist")
    c1.Update()

    c2 = ROOT.TCanvas("c2","Subassembly 2",1200,700)
    c2.Divide(12,6)
    for pmt in range(72,144):
        c2.cd(pmt+1-72)
        histos[pmt].Draw("hist")
    c2.Update()

    c3 = ROOT.TCanvas("c3","Subassembly 3",1200,700)
    c3.Divide(12,6)
    for pmt in range(144,216):
        c3.cd(pmt+1-144)
        histos[pmt].Draw("hist")
    c3.Update()

    c4 = ROOT.TCanvas("c4","Subassembly 4",1200,700)
    c4.Divide(12,6)
    for pmt in range(216,288):
        c4.cd(pmt+1-216)
        histos[pmt].Draw("hist")
    c4.Update()

    input("Please press enter to close images.")  

#plot_pmt(0,0)
#plot_event(11)
#plot_event(1)

print(dig_adc.shape)
print(dig_adc[0][0][0])

save_path = '/home/skbarcus/JLab/SBS/HCal/Machine_Learning/G4SBS_digitization/data_arrays/'
print('Saving numpy arrays to',save_path,'.')
if save_data == 1:
    np.save(save_path+"fnucl_"+rootfile[:-5], fnucl)
    np.save(save_path+"adc_int_"+rootfile[:-5], adc_int)
    np.save(save_path+"pmt_hits_"+rootfile[:-5], pmt_hits)

print("The script took %.2f seconds (%.2f minutes or %.2f hours) to run." % (time.time() - start, (time.time() - start)/60.,(time.time() - start)/60./60.)) #Print time to run.
