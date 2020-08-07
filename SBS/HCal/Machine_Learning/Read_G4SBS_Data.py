from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
from mpl_toolkits.axes_grid1.mpl_axes import Axes
import numpy as np
import ROOT
#from ROOT import TChain, TSelector, TTree
import sys

from keras.utils import np_utils #Utilities to transform data.
from sklearn.preprocessing import normalize
#from root_numpy import root2array, testdata
#root2array(testdata.get_filepath('single1.root'))[:20]

import time
start = time.time() #Start a timer.

save_data = 0
read_all = 1
loop = 30000
test_evt = 2#31823
nsamps = 20

data = "/home/skbarcus/JLab/SBS/HCal/Analysis/Simulation/G4SBS/rootfiles/"
#infile = sys.argv[1] #Use file specified when code is run. 

chain = ROOT.TChain("T")
for i in range(0,20):
    rootfile = "gmn_kin07_r"+("%d" % i)+".root"
    chain.AddFile(data+rootfile)
    print("Adding rootfile ",data+rootfile," to chain.")

print(chain.GetEntries())
if read_all == 1:
    loop = chain.GetEntries()

fnucl = []         #0 = neutron, 1 = proton.
front_nhits = []   #Hits on the front panel of HCal.
nhits = []         #Scintillator hits.
pmt = []           #PMT module with a hit.
pmt_edep = []      #Energy deposited in a single PMT (scintillator).
tavg = []          #Average time scint hit occurred.

#Loop over entries.
for entryNum in range (0, loop):
    chain.GetEntry(entryNum)
    fnucl_evt = getattr(chain,"fnucl")#0 = neutron, 1 = proton.
    fnucl.append(fnucl_evt)
    front_nhits_evt = getattr(chain,"Harm.HCalFrontPlate.hit.nhits")
    front_nhits.append(front_nhits_evt)
    nhits_evt = getattr(chain,"Harm.HCalScint.hit.nhits")
    nhits.append(nhits_evt)
    pmt_evt = getattr(chain,"Harm.HCalScint.hit.cell")
    pmt_evt = np.array(pmt_evt)
    #print(pmt_evt)
    pmt.append(pmt_evt)
    pmt_edep_evt = getattr(chain,"Harm.HCalScint.hit.sumedep")
    pmt_edep_evt = np.array(pmt_edep_evt)
    pmt_edep.append(pmt_edep_evt)
    tavg_evt = getattr(chain,"Harm.HCalScint.hit.tavg")
    tavg_evt = np.array(tavg_evt)
    tavg.append(tavg_evt)
    #fnucl_evt = getattr(chain,"Harm.HCalScint.hit.yhit")
    if entryNum%10000==0:
        print("Read in %d events." % entryNum)
    
fnucl = np.array(fnucl)
print(fnucl)
front_nhits = np.array(front_nhits)
print(front_nhits)
nhits = np.array(nhits)
print(nhits)
pmt = np.array(pmt)
print(pmt, " PMT data type = ", type(pmt))
pmt_edep = np.array(pmt_edep)
print(pmt_edep)
tavg = np.array(tavg)
print(tavg)
#print(pmt[19][22])

if save_data == 1:
    np.save(data + "g4sbs_kin7_fnucl", fnucl)
    np.save(data + "g4sbs_kin7_front_nhits", front_nhits)
    np.save(data + "g4sbs_kin7_nhits",nhits )
    np.save(data + "g4sbs_kin7_pmt", pmt)
    np.save(data + "g4sbs_kin7_pmt_edep", pmt_edep)
    np.save(data + "g4sbs_kin7_tavg", tavg)

print("The script took %.2f seconds (%.2f minutes or %.2f hours) to run." % (time.time() - start, (time.time() - start)/60.,(time.time() - start)/60./60.)) #Print time to run.
