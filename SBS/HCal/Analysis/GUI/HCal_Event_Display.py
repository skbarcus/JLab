from tkinter import Tk, Label, Button, StringVar, W, E, Frame, RIGHT, BOTH, RAISED, LEFT, CENTER, TOP, BOTTOM, Canvas, Scrollbar, Radiobutton
from tkinter.ttk import Style, Entry
import json
import os
from datetime import date, datetime

import Lib_HCal_Event_Display as fn
import ROOT
import sys
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
from mpl_toolkits.axes_grid1.mpl_axes import Axes

read_all = 0
loop = 1
test_evt = 2#31823

data = "/home/skbarcus/Machine_Learning/Tutorials/Data/"

#infile = sys.argv[1]
#infile = "/home/skbarcus/JLab/SBS/HCal/Analysis/rootfiles/fadc_f1tdc_1710.root"
infile = "/home/skbarcus/JLab/SBS/HCal/Analysis/Cosmics/rootfiles/fadc_f1tdc_820.root"

print("Reading from", infile)

inFile = ROOT.TFile.Open(infile," READ ")

tree = inFile.Get("T")

h1 = ROOT.TH1D("h1","h1 test",20,0,50000)

print(tree.GetEntries())

#for entryNum in range (0, loop):
    #tree.GetEntry(entryNum)

tree.GetEntry(0)
nsamps = getattr(tree,"sbs.hcal.nsamps")
print('The number of fADC samples is '+str(nsamps[0])+'.')

if read_all == 1:
    loop = tree.GetEntries()

adc_samps_vals = []
adc_vals = []
tdc_vals = []
hit_vals = []

class DVCS_Pulser_Control_GUI:
    def __init__(self, primary):
        #super().__init__(primary)
        self.primary = primary
        primary.geometry("800x700")
        primary.title("DVCS Pulser Control GUI")

        #Create a frame to hold the exit button.
        exit_frame = Frame(primary, relief=RAISED, borderwidth=2)
        exit_frame.pack(side='bottom')
        exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=root.quit, bg = "red")
        exit_btn.pack(side="bottom")

        #Create a frame to hold directions and settings along with the save settings button.
        info_frame = Frame(primary, relief=RAISED, borderwidth=2)
        info_frame.pack(side='left', fill=BOTH, expand=True)

        #Create a label with the instructions to use the control GUI.
        instructions_label_text = StringVar()
        instructions_label_text.set("This GUI controls the DVCS pulser which\nis used to pulse the LEDs in HCal\nfor gain monitoring purposes. This pulser is\nused in conjunction with a NIM gate generator\nto produce logical pulses compatible with\nthe HCal LED power boxes.\n\n")
        instructions_label = Label(info_frame, width = 40,textvariable=instructions_label_text,font='Helvetica 12 bold')
        instructions_label.pack(side="top")

        #Create frame to set the individual channel voltage outputs.
        display_frame = Frame(primary, relief=RAISED, borderwidth=2)
        display_frame.pack(side='left', fill=BOTH, expand=True)

        #Create an entry for the event number to be displayed.
        display_entry = Entry(info_frame,width=8)
        display_btn = Button(info_frame, text='Display\nEvent', width=8, height=4, font='Helvetica 8 bold', command=lambda info_frame=info_frame, display_frame=display_frame, tree=tree, display_entry=display_entry: fn.display_event(info_frame,display_frame,tree,display_entry), bg = "grey")
        display_btn.pack(side="bottom")

        display_entry.pack(side="bottom")
        display_entry.insert(0,0)

        def key(event):
            print("pressed", repr(event.char))

        def callback(event):
            print("clicked at", event.x, event.y)

        canvas= Canvas(display_frame, width=100, height=100, bd=2, relief=RAISED)
        canvas.bind("<Key>", key)
        canvas.bind("<Button-1>", callback)
        canvas.pack()

root = Tk()
my_gui = DVCS_Pulser_Control_GUI(root)
root.mainloop()
