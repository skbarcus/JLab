from tkinter import Tk, Label, Button, StringVar, W, E, Frame, RIGHT, BOTH, RAISED, LEFT, CENTER, TOP, BOTTOM, Canvas, Scrollbar
from tkinter.ttk import Style, Entry
import json

nrows = 24
ncols = 12
channels = nrows * ncols
pmt_mods = []    #Holds PMT module numbers.
buttons = []

#Dictionary containing all of the relevant connections for each PMT.
# reading the data from the file 
with open('hcal_connections.json') as f: 
    connections = f.read()
#print(type(connections))
# reconstructing the data as a dictionary 
connections = json.loads(connections)
#print(type(connections))
#print(connections['1'][0])

#Fill PMT module numbers.
for i in range(0,channels):
    pmt_mods.append(i+1)
    #print(pmt_mods[i])

#Fill a list with the names for the buttons.
button_names = []
for i in range(0,channels):
    button_names.append("Button"+str(i+1))
#print(button_names)

class MyFirstGUI:
    LABEL_TEXT = [
        "This is our first GUI!",
        "Actually, this is our second GUI.",
        "We made it more interesting...",
        "...by making this label interactive.",
        "Go on, click on it again.",
    ]
    def __init__(self, primary):
        self.primary = primary
        primary.title("Hadron Calorimeter")

        Style().configure("TButton", padding=(0, 5, 0, 5),font='serif 10')

        #frame1 = Frame(primary)

        #exit_frame = Frame(primary, relief=RAISED, borderwidth=2)
        #exit_frame.pack(side='bottom')
        #exit_btn = Button(exit_frame, text='Close', width=2, height=1, font='Helvetica 8 bold', command=root.destroy, bg = "red")
        #exit_btn.pack(side="bottom")

        for i in range(0,ncols):
            primary.columnconfigure(i, pad=3)

        for i in range(0,nrows+1):
            primary.rowconfigure(i, pad=3)

        #Fill a list with the actual buttons.
        #buttons = []
        for i in range(0,channels):
            self.button = Button(primary, text=pmt_mods[i], width=1, height=1)#, command=self.fconnections)
            #button["command"] = self.fconnections
            self.button.bind("<Button-1>", self.fconnections)
            buttons.append(self.button)

        for i in range(0,nrows):
            for j in range(0,ncols):
                buttons[i*ncols+j].grid(row=i, column=j)
                #self.buttons[i*ncols+j].pack()

        #exit_btn = Button(primary, text='Close', width=2, height=1, font='Helvetica 8 bold', command=root.destroy, bg = "red") #Don't use destroy! It hangs the terminal for a long time!
        exit_btn = Button(primary, text='Close', width=2, height=1, font='Helvetica 8 bold', command=root.quit, bg = "red")
        exit_btn.grid(row=24, column=5, columnspan=2)

        #test_button = Button(primary, text='Test', command=self.fconnections)
        #self.test_button.pack()

        #self.label_index = 0
        #self.label_text = StringVar()
        #self.label_text.set(self.LABEL_TEXT[self.label_index])
        #self.label = Label(primary, textvariable=self.label_text)
        #self.label.bind("<Button-1>", self.cycle_label_text)
        #self.label.pack()

        #self.greet_button = Button(primary, text="Greet", command=self.greet)
        #self.greet_button.pack()

        #self.close_button = Button(primary, text="Close", command=primary.quit)
        #self.close_button.pack()

    def greet(self):
        print("Greetings!")

    def cycle_label_text(self, event):
        print(self,'  ',event)
        print(event.widget)
        self.label_index += 1
        self.label_index %= len(self.LABEL_TEXT) # wrap around
        self.label_text.set(self.LABEL_TEXT[self.label_index])

    def old_fconnections(self, event):
        pmt = event.widget.cget('text')
        print('PMT Module ',pmt,': amplifier channel = ', connections[pmt][0],', fADC channel = ', connections[pmt][1], ', HV channel = ',connections[pmt][2])
        #print(event.widget.cget('text'))

    def fconnections(self, event):
        pmt = event.widget.cget('text')
        pmt = str(pmt)
        print('******************** Information for PMT Module '+pmt+' ********************')
        print('PMT module '+pmt+' is powered by HV channel '+str(connections[pmt][11])+'.')
        print('PMT module '+pmt+'\'s output goes to amplifier '+str(connections[pmt][0])+'.')
        print('This signal terminates at fADC '+str(connections[pmt][3])+', F1TDC '+str(connections[pmt][9])+', and summing module '+str(connections[pmt][10])+'.')
        print('The fADC data flow follows: amplfier '+str(connections[pmt][0])+' --> front-end fADC patch panel '+str(connections[pmt][1])+' --> DAQ fADC patch panel'+str(connections[pmt][2])+' --> fADC '+str(connections[pmt][3])+'.')
        print('The TDC data flow follows: amplfier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> front-end f1TDC discriminator '+str(connections[pmt][5])+' --> front-end TDC patch panel '+str(connections[pmt][6])+' --> DAQ TDC patch panel '+str(connections[pmt][7])+' --> DAQ TDC discriminator '+str(connections[pmt][8])+' --> F1TDC '+str(connections[pmt][9])+'.')
        print('The summing module data flow follows: amplifier '+str(connections[pmt][0])+' --> splitter panel '+str(connections[pmt][4])+' --> summing module '+str(connections[pmt][10])+'.')
        #print(event.widget.cget('text'))

root = Tk()
#print(root.call("info", "patchlevel")) #Prints tkinter version.
my_gui = MyFirstGUI(root)
root.mainloop()
