from tkinter import Tk, Label, Button, StringVar, W, E, Frame, RIGHT, BOTH, RAISED, LEFT, CENTER, TOP, BOTTOM, Canvas, Scrollbar
from tkinter.ttk import Style, Entry

class DAC_Control_GUI:
    def __init__(self, primary):
        #super().__init__(primary)
        self.primary = primary
        primary.geometry("800x600")
        primary.title("VME Digital to Analog Converter Control GUI")

        #Create a frame to hold the exit button.
        exit_frame = Frame(primary, relief=RAISED, borderwidth=2)
        exit_frame.pack(side='bottom')
        exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=root.quit, bg = "red")
        exit_btn.pack(side="bottom")

        #Create a frame to hold directions and settings along with the save settings button.
        info_frame = Frame(primary, relief=RAISED, borderwidth=2)
        info_frame.pack(side='left', fill=BOTH, expand=True)

        #self.test_name = StringVar()
        self.test_entry = Entry(info_frame, width = 5)
        #self.test_entry = Entry(info_frame, width = 5, textvariable = self.test_name)
        #test_entry.insert(0,'0.00')
        self.test_entry.pack(side="top")

        save_btn = Button(info_frame, text='Save Settings', width=10, height=4, font='Helvetica 12 bold', command=self.greet, bg = "yellow")
        save_btn.pack(side="bottom")

        #Create frame to set the individual channel voltage outputs.
        ch_settings_frame = Frame(primary, relief=RAISED, borderwidth=2)
        ch_settings_frame.pack(side='left', fill=BOTH, expand=True)

        nchs = 16                  #Number of VME DAC output channels.
        nrows = nchs
        ncols = 4
        self.ch_labels = []             #Array to hold channel labels.
        self.readback_labels = []      #Array to hold current voltage setting readback labels.
        self.predicted_threshold_labels = []   #Array to hold predicted VME DAC threshold setting.
        self.voltage_setting_entries = []      #Array to hold entry widgets for entering the DAC voltage settings.

        entry_names = []

        for ch in range(nchs):
            #Create and store the channel labels.
            self.ch_label_text = StringVar()
            self.ch_label_text.set("Channel "+str(ch))
            self.ch_label = Label(ch_settings_frame, width = 9,textvariable=self.ch_label_text)
            self.ch_labels.append(self.ch_label)

            #Create and store the channel voltage readback labels.
            self.readback_label_text = StringVar()
            self.readback_label_text.set("0.00")
            self.readback_label = Label(ch_settings_frame, width = 9,textvariable=self.readback_label_text)
            self.readback_labels.append(self.readback_label)

            #Create and store the predicted discriminator threshold labels.
            self.thr_label_text = StringVar()
            self.thr_label_text.set("0.00")
            self.thr_label = Label(ch_settings_frame, width = 9,textvariable=self.thr_label_text)
            self.predicted_threshold_labels.append(self.thr_label)

            #Create and store the entry widgets for setting the channel voltage outputs.
            self.voltage_setting_entry = Entry(ch_settings_frame, width = 5)
            self.voltage_setting_entry.insert(0,'0.00')
            self.voltage_setting_entries.append(self.voltage_setting_entry)

        #print(voltage_setting_entries[1])
        #Create labels for the columns.
        col_labels = []
        col_labels_text = ['Output Channels','Current Voltage\n Setting (Volts)','Predicted\n Discriminator\n Threshold (mV)','Enter New\n Voltage Setting\n (Volts)']
        for col in range(ncols):
            col_label_text = StringVar()
            col_label_text.set(col_labels_text[col])
            col_label = Label(ch_settings_frame, width = 9,textvariable=col_label_text,font='Helvetica 10 bold')
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
            for col in range(0,ncols):
                if row==0:
                    col_labels[col].grid(row=row, column=col, sticky='NSEW')

    def greet(self):
        print("Greetings!")
        print(self.voltage_setting_entries[0].get())
        print(self.test_entry.get())

"""
        var = StringVar()
        var.set("Channel 0")
        ch0_label = Label(ch_settings_frame, width = 9,textvariable=var)
        ch0_label.pack(side='left',padx = 5, pady = 5)

        ch0_set = Entry(ch_settings_frame, width = 5)
        ch0_set.insert(0,'0.0')
        ch0_set.pack(side='left',padx = 5, pady = 5)
"""

root = Tk()
#ex = Example()
#root.geometry("800x600")
my_gui = DAC_Control_GUI(root)
root.mainloop()
