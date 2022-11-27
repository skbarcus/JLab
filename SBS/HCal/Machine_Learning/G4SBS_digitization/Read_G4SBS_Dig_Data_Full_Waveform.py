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

save_data = 0 #Careful! 50k events with nsamps=20 -> 2.1 G numpy array for full ADC.
read_all = 0  #Careful! 50k events with nsamps=20 -> 17.5 min to run code.
max_evts = 15
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
else:
    loop = max_evts

print(rootfile,"has ",loop,"entries.")

fnucl = []         #0 = neutron, 1 = proton.
front_nhits = []   #Hits on the front panel of HCal.
nhits = []         #Scintillator hits.
pmt = []           #PMT module with a hit.
pmt_edep = []      #Energy deposited in a single PMT (scintillator).
tavg = []          #Average time scint hit occurred.
dig_adc = []       #Holds digitized fADC values in an array of size dig_adc[nsamps] with dig_adc[0] holding all ADC values for the first fADC bin in an event.
new_dig_adc = []   #Holds digitized fADC values in an array of size dig_adc[nsamps] with dig_adc[0] holding all ADC values for the first fADC bin in an event. No crazy shape this time and full nPMT entries.
pmt_adc = []       #Holds the fADC values for each PMT in each event. If G4SBS said PMT didn't fire, it fills the PMT with random pedestal noise.
pmt_hits = []      #1 if the PMT fired according to G4SBS and 0 if it didn't.
adc_int = []       #Array to hold all 288 fADC values for each event.
new_adc_int = []   #Array to hold all 288 fADC values for each event built from new_dig_adc.
tdc = []
for i in range(0,nsamps):
    dig_adc.append([])
srau = []          #Integrals of PMT fADCs.
chan = []          #Stores the PMT channel of the hit. Values 0-287.

#Initialized digitized fADC values array with zeros. new_dig_adc[event][pmt][samp].
new_dig_adc = np.zeros((loop,npmts,nsamps))

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
#Loop over all events.
for event in range(0,loop):
    pmt_hits_evt = []
    fadc_integrals = []
    hit_idx = 0 #Index of the PMT in the dig_adc array. i.e. [pmt 12, pmt 13, pmt 14]->pmt 13 hit_idx=1.

    #Loop over all PMTs. 
    for pmt in range(0,npmts):
        fadc_integral = 0
        if pmt in chan[event]:
            pmt_hits_evt.append(1)#Fill PMT hits.
            #Loop over all fADC samples.
            for samp in range(0,nsamps):
                #Fill the fADC samples with the values from simulation.
                new_dig_adc[event][pmt][samp]=dig_adc[samp][event][hit_idx]
                #Sum the fADC integral for this PMT.
                fadc_integral = fadc_integral+new_dig_adc[event][pmt][samp]
            fadc_integrals.append(fadc_integral)
            hit_idx=hit_idx+1
        else:
            pmt_hits_evt.append(0)#Fill PMT misses.
            #Loop over all fADC samples.
            for samp in range(0,nsamps):
                #Fill the fADC samples with pedestal noise for PMTs with no hit.
                new_dig_adc[event][pmt][samp]=round(gauss(300, 1.6))
                #Sum the fADC integral for this PMT.
                fadc_integral = fadc_integral+new_dig_adc[event][pmt][samp]
            fadc_integrals.append(fadc_integral)

    #Append the 288 PMT hit/miss values to the hits array.
    pmt_hits.append(pmt_hits_evt)
    #Append the fADC integral values.
    new_adc_int.append(fadc_integrals)

    #Print some timing and progress info.
    if event%5000==0:
        print('Building data arrays is '+str((event/loop)*100.)+'% complete.')
    elif event==(loop-1):
        print('Building data arrays is 100% complete.')
        print("This took %.2f seconds (%.2f minutes or %.2f hours)." % (time.time() - start, (time.time() - start)/60.,(time.time() - start)/60./60.))

pmt_hits = np.array(pmt_hits)
#print('pmt_hits',pmt_hits)

new_dig_adc = np.array(new_dig_adc)
print('type(new_dig_adc)',type(new_dig_adc))
print('new_dig_adc.shape',new_dig_adc.shape)
print('new_dig_adc[0]',new_dig_adc[0][228])

new_adc_int = np.array(new_adc_int)
print('type(new_adc_int)',type(new_adc_int))
print('new_adc_int.shape',new_adc_int.shape)
print('new_adc_int[0]',new_adc_int[0][228])

#print(len(fnucl))
#print(fnucl)

def plot_pmt(event,pmt):
    xaxis = []
    for i in range(0,nsamps):
        xaxis.append(i)
    print("xaxis =",xaxis)

    yaxis = []
    for samp in range(0,nsamps):
        yaxis.append(new_dig_adc[event][pmt][samp])
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
    #Create and fill root histos for each PMT's fadc samples for a single hit event. Maybe add cuts on eng later.
    #for hit in range(0,len(dig_adc[0][event])):
    for pmt in range(0,npmts):
        #pmt = chan[event][hit]
        histos.append(ROOT.TH1F('hfadc%d' % (pmt+1),'PMT %d' % (pmt+1),nsamps,0,nsamps))
        #histos[pmt].SetTitleSize(10.,"t");
        ROOT.gStyle.SetTitleSize(1,"t");
        ROOT.gStyle.SetTitleY(1.5)
        histos[pmt].SetStats(0);
        histos[pmt].GetYaxis().SetLabelSize(0.15);
        histos[pmt].GetYaxis().SetLabelOffset(-0.2);
        histos[pmt].GetYaxis().SetNdivisions(5);
        histos[pmt].SetLineWidth(2)
        xaxis = []
        yaxis = []
        #Color the PMTs with hits red.
        if pmt in chan[event]:
            histos[pmt].SetLineColor(2)
        for samp in range(0,nsamps):
            xaxis.append(samp)
            yaxis.append(new_dig_adc[event][pmt][samp])
            histos[pmt].Fill(xaxis[samp],yaxis[samp])
            #print("xaxis =",xaxis)
            #print("yaxis =",yaxis)

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

#plot_pmt(0,240)
#plot_event(11)
plot_event(0)

print(new_dig_adc.shape)
print(new_dig_adc[0][228][0])

save_path = '/home/skbarcus/JLab/SBS/HCal/Machine_Learning/G4SBS_digitization/data_arrays/'
print('Saving numpy arrays to',save_path,'.')
if save_data == 1:
    np.save(save_path+"fnucl_"+rootfile[:-5], fnucl)
    np.save(save_path+"adc_int_"+rootfile[:-5], new_adc_int)
    np.save(save_path+"adc_full_"+rootfile[:-5], new_dig_adc)
    np.save(save_path+"pmt_hits_"+rootfile[:-5], pmt_hits)

print("The script took %.2f seconds (%.2f minutes or %.2f hours) to run." % (time.time() - start, (time.time() - start)/60.,(time.time() - start)/60./60.)) #Print time to run.
