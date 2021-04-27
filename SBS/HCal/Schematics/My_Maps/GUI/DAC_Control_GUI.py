from tkinter import Tk, Label, Button, StringVar, W, E, Frame, RIGHT, BOTH, RAISED, LEFT, CENTER, TOP, BOTTOM, Canvas, Scrollbar
from tkinter.ttk import Style, Entry
import json
import os

class DAC_Control_GUI:
    def __init__(self, primary):
        #super().__init__(primary)
        self.primary = primary
        primary.geometry("800x700")
        primary.title("VME Digital to Analog Converter Control GUI")

        global save_file
        save_file = 'mpv904_DAC_Voltage_Settings.json'

        #os.system('ls')#Reminder for how to run commands in terminal with Python.

        #Read in json file containing a list with all of the saved voltage channel settings.
        voltage_settings = get_voltage_settings(save_file)
        for ch in range(16):
            print('Output channel',ch,'is set to',voltage_settings[ch],'Volts.')

        #Create a frame to hold the exit button.
        exit_frame = Frame(primary, relief=RAISED, borderwidth=2)
        exit_frame.pack(side='bottom')
        exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=root.quit, bg = "red")
        exit_btn.pack(side="bottom")

        #Create a frame to hold directions and settings along with the save settings button.
        info_frame = Frame(primary, relief=RAISED, borderwidth=2)
        info_frame.pack(side='left', fill=BOTH, expand=True)

        save_btn = Button(info_frame, text='Save Settings', width=10, height=4, font='Helvetica 12 bold', command=self.save_settings, bg = "yellow")
        save_btn.pack(side="bottom")

        #canvas = Canvas(info_frame, width = 60,height = 60)
        #my_oval = canvas.create_oval(5, 5, 30, 30)#(x0,y0,x1,y1)=(top left corner of rectangle, bottom right corner)
        #canvas.itemconfig(my_oval, fill="red")
        #canvas.create_text(17,18,font='Helvetica 8 bold',text='OFF')
        #canvas.pack(side='top')

        #Create frame to set the individual channel voltage outputs.
        ch_settings_frame = Frame(primary, relief=RAISED, borderwidth=2)
        ch_settings_frame.pack(side='left', fill=BOTH, expand=True)

        self.nchs = 16                         #Number of VME DAC output channels.
        nrows = self.nchs                      #Number of rows of channels settings.
        ncols = 6                              #Number of columns per channel setting row.
        self.ch_labels = []                    #Array to hold channel labels.
        self.readback_labels = []              #Array to hold current voltage setting readback labels.
        self.predicted_threshold_labels = []   #Array to hold predicted VME DAC threshold setting.
        self.voltage_setting_entries = []      #Array to hold entry widgets for entering the DAC voltage settings.
        self.power_buttons = []                #Array to hold the power on/off buttons for each output channel.
        self.indicator_lights = []             #Array to hold the indicator lights to see if a channel is energized.

        entry_names = []

        for ch in range(self.nchs):
            #Create and store the channel labels.
            self.ch_label_text = StringVar()
            self.ch_label_text.set("Channel "+str(ch))
            self.ch_label = Label(ch_settings_frame, width = 9,textvariable=self.ch_label_text)
            self.ch_labels.append(self.ch_label)

            #Create and store the channel voltage readback labels.
            self.readback_label_text = StringVar()
            self.readback_label_text.set(str(voltage_settings[ch]))
            self.readback_label = Label(ch_settings_frame, width = 9,textvariable=self.readback_label_text)
            self.readback_labels.append(self.readback_label)

            #Create and store the predicted discriminator threshold labels.
            self.thr_label_text = StringVar()
            self.thr_label_text.set("NA")
            self.thr_label = Label(ch_settings_frame, width = 9,textvariable=self.thr_label_text)
            self.predicted_threshold_labels.append(self.thr_label)

            #Create and store the entry widgets for setting the channel voltage outputs.
            self.voltage_setting_entry = Entry(ch_settings_frame, width = 5)
            self.voltage_setting_entry.insert(0,str(voltage_settings[ch]))
            self.voltage_setting_entries.append(self.voltage_setting_entry)

            #Create and store the buttons for enabling/disabling the channel voltage outputs.
            self.power_btn = Button(ch_settings_frame, text='Enable\n Voltage', width=3, height=1, font='Helvetica 8 bold')
            self.power_btn.configure(command=lambda button=self.power_btn, buttons=self.power_buttons, indicator_lights=self.indicator_lights: toggle_power(button,buttons,indicator_lights))
            self.power_buttons.append(self.power_btn)

            #Create and store indicator 'lights' to show whether a channel is energized or not.
            self.canvas = Canvas(ch_settings_frame, width = 45,height = 45)
            my_oval = self.canvas.create_oval(15, 5, 40, 30)#(x0,y0,x1,y1)=(top left corner of rectangle, bottom right corner)
            self.canvas.itemconfig(my_oval, fill="red")
            self.canvas.create_text(27,18,font='Helvetica 8 bold',text='OFF')
            self.indicator_lights.append(self.canvas)

        #Create labels for the columns.
        col_labels = []
        col_labels_text = ['Output Channels','Current Voltage\n Setting (Volts)','Predicted\n Discriminator\n Threshold (mV)','Enter New\n Voltage Setting\n (Volts)','Enable/Disable\n Channel Voltage','Power\n Indicator']
        for col in range(ncols):
            col_label_text = StringVar()
            col_label_text.set(col_labels_text[col])
            col_label = Label(ch_settings_frame, width = 14,height=3,textvariable=col_label_text,font='Helvetica 10 bold')
            col_labels.append(col_label)
            

        #Define the rows and cols for the button grid.
        for row in range(0,nrows+1):
            ch_settings_frame.rowconfigure(row, pad=3, weight=1)
            
        for col in range(0,ncols):
            ch_settings_frame.columnconfigure(col, pad=3, weight=1)

        #Place widgets on the grid.
        for row in range(0,nrows+1):
            if row!=0:
                self.ch_labels[row-1].grid(row=row, column=0, sticky='NSEW')
                self.readback_labels[row-1].grid(row=row, column=1, sticky='NSEW')
                self.predicted_threshold_labels[row-1].grid(row=row, column=2, sticky='NSEW')
                self.voltage_setting_entries[row-1].grid(row=row, column=3)
                self.power_buttons[row-1].grid(row=row, column=4)
                self.indicator_lights[row-1].grid(row=row, column=5)
            for col in range(0,ncols):
                if row==0:
                    col_labels[col].grid(row=row, column=col, sticky='NSEW')

    def save_settings(self):
        voltage_settings = [] #Holds the voltage settings for each of the 16 channels and reads them from the entry boxes when the save button is clicked.
        save_file = 'mpv904_DAC_Voltage_Settings.json'
        for ch in range(self.nchs):
            voltage_settings.append(self.voltage_setting_entries[ch].get())
            print('Output channel',ch,'is now set to',voltage_settings[ch],'Volts.')
        print('New voltage settings saved to',save_file+'.')
        #Write dictionary to json file.
        out_file = json.dumps(voltage_settings)
        f = open(save_file,"w")
        f.write(out_file)
        f.close()

#Function to return voltage settings.
def get_voltage_settings(save_file):
    #Read in json file containing a list with all of the saved voltage channel settings.
    with open(save_file) as f: 
        voltage_settings = f.read()
        voltage_settings = json.loads(voltage_settings)
    return voltage_settings

#Function to enable and disable the voltage output for individual channels.
def toggle_power(button,buttons,indicator_lights):
    voltage_settings = get_voltage_settings(save_file)
    text = button.cget('text')
    off = 'Enable\n Voltage'
    on = 'Disable\n Voltage'
    
    #Check which button was pressed by searching the array holding the power buttons. 
    if button in buttons:
        ch = buttons.index(button)
        voltage = voltage_settings[ch]
        print('Channel',ch,'is set to',voltage,'Volts.')

    if text==off:
        print('Turning channel',ch,'on.')
        button['text']=on
        my_oval = indicator_lights[ch].create_oval(15, 5, 40, 30)#(x0,y0,x1,y1)=(top left corner of rectangle, bottom right corner)
        indicator_lights[ch].itemconfig(my_oval, fill="green")
        indicator_lights[ch].create_text(27,18,font='Helvetica 8 bold',text='ON')

    if text==on:
        print('Turning channel',ch,'off.')
        button['text']=off
        my_oval = indicator_lights[ch].create_oval(15, 5, 40, 30)#(x0,y0,x1,y1)=(top left corner of rectangle, bottom right corner)
        indicator_lights[ch].itemconfig(my_oval, fill="red")
        indicator_lights[ch].create_text(27,18,font='Helvetica 8 bold',text='OFF')

root = Tk()
my_gui = DAC_Control_GUI(root)
root.mainloop()