#Library of functions for HCal_Event_Display.py.
import ROOT
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
from mpl_toolkits.axes_grid1.mpl_axes import Axes
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)

def test():
    print('Test')

def test_command():
    print('Test.')

def display_event(info_frame,display_frame,tree,display_entry):
    event = display_entry.get()
    print('Now displaying event '+str(event)+'.')

    event = int(event)
    tree.GetEntry(event)
    nsamps = getattr(tree,"sbs.hcal.nsamps")
    nsamps = int(nsamps[0])
    samps = getattr(tree,"sbs.hcal.samps")
    print('The number of fADC samples is '+str(nsamps)+'.')

    xaxis = []
    xaxis.clear()
    for i in range(nsamps):
        xaxis.append(i)

    yaxis = []
    yaxis.clear()
    for i in range(nsamps):
        yaxis.append(samps[i])

    #Plot 1D ROOT histogram of a single PMT fADC channel.
    #h1dfadc = ROOT.TH1F("hfadc","hfadc test",20,0,20)
    #for i in range (0,nsamps):
        #h1dfadc.Fill(xaxis[i],yaxis[i])

    #h1dfadc.Draw("hist")
    #input("Please press enter to close images.") 

    fig, ax = plt.subplots(figsize=(4,4))
    hist1 = ax.hist(xaxis,bins=nsamps,weights=yaxis,histtype='bar',ec='black')
    ax.set_ylim(200,500)
    ax.set_xlim(0,nsamps-1)
    #plt.show()

    #canvas.get_tk_widget().destroy()

    #try: 
        #canvas.get_tk_widget().destroy()
    #except:
        #print('Error!')
        #pass 

    # creating the Tkinter canvas
    # containing the Matplotlib figure
    print(display_frame.winfo_children())
    canvas = FigureCanvasTkAgg(fig, master=display_frame)  

    #try: 
        #canvas.get_tk_widget().pack_forget()
    #except AttributeError: 
        #pass 
    #canvas.delete("all")

    #canvas.get_tk_widget().destroy()

    canvas.draw()

    #canvas.get_tk_widget().destroy()

    # placing the canvas on the Tkinter window
    canvas.get_tk_widget().pack()

    #canvas.get_tk_widget().destroy()
    #canvas.get_tk_widget().forget_pack()

    # creating the Matplotlib toolbar
    toolbar = NavigationToolbar2Tk(canvas,display_frame)
    toolbar.update()
    
    # placing the toolbar on the Tkinter window
    canvas.get_tk_widget().pack()

    #canvas.get_tk_widget().destroy()

    #ax.clear()
    #canvas.get_tk_widget().forget_pack()
    #canvas.get_tk_widget().delete("all")
    #canvas.delete("all")
