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
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, 
NavigationToolbar2Tk)
from random import seed
from random import random

test_var = 0

class Canvas_Destroy_Test:
    def __init__(self, primary):
        #super().__init__(primary)
        self.primary = primary
        primary.geometry("800x700")
        primary.title("Canvas_Destroy_Test")

        xaxis = []
        xaxis.clear()
        for i in range(10):
            xaxis.append(i)

        yaxis = []
        yaxis.clear()
        for i in range(10):
            yaxis.append(random())

        fig, ax = plt.subplots(figsize=(4,4))
        hist1 = ax.hist(xaxis,bins=10,weights=yaxis,histtype='bar',ec='black')
        #ax.set_ylim(200,500)
        #ax.set_xlim(0,nsamps-1)
        #plt.show()
        
        canvas = FigureCanvasTkAgg(fig, master=root)
        #canvas.get_tk_widget().bind("<Key>", key)
        #canvas.get_tk_widget().bind("<Button-1>", callback)
        canvas.get_tk_widget().bind("<Button-1>", popupBonus)
        #canvas.bind("<Key>", key)
        #canvas.bind("<Button-1>", callback)
        #canvas.draw()
        #canvas.get_tk_widget().pack()
        
        #toolbar = NavigationToolbar2Tk(canvas,root)
        #toolbar.update()
        #canvas.get_tk_widget().pack()

        btn = Button(root, text='Draw Plot', width=8, height=4, font='Helvetica 8 bold', command=lambda canvas=canvas, ax=ax: plot(canvas,ax), bg = "grey")
        #btn.pack(side="bottom")
        btn.grid(row=0,column=0)

def key(event):
    print("pressed", repr(event.char))

def callback(event):
    print("clicked at", event.x, event.y)

def popupBonus(event):
    popupBonusWindow = Tk()
    popupBonusWindow.wm_title("Window")
    labelBonus = Label(popupBonusWindow, text="Input")
    labelBonus.grid(row=0, column=0)
    B1 = Button(popupBonusWindow, text="Close", width=8, height=4, command=popupBonusWindow.destroy)
    B1.grid(row=0,column=0)

def plot(canvas,ax):
    #try: 
    #    canvas.get_tk_widget().destroy()
    #except:
    #    print('Error!')
    #    pass 

    print('Plotting')
    print(root.winfo_children())

    xaxis = []
    xaxis.clear()
    for i in range(10):
        xaxis.append(i)

    yaxis = []
    yaxis.clear()
    for i in range(10):
        yaxis.append(random())

    #fig, ax = plt.subplots(figsize=(4,4))
    ax.clear()
    hist1 = ax.hist(xaxis,bins=10,weights=yaxis,histtype='bar',ec='black')

    canvas.draw()
    #canvas.get_tk_widget().pack()
    canvas.get_tk_widget().grid(row=0,column=1)

    global test_var
    print(test_var)
    test_var = 1
    print(test_var)

root = Tk()
my_gui = Canvas_Destroy_Test(root)
root.mainloop()
