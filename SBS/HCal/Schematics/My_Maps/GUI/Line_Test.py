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

class Example(Frame):

    def __init__(self):
        super().__init__()

        self.initUI()


    def initUI(self):

        self.master.title("Lines")
        self.pack(fill=BOTH, expand=1)

        canvas = Canvas(self)
        canvas.create_line(15, 25, 200, 25)
        canvas.create_line(300, 35, 300, 200, dash=(4, 2))
        canvas.create_line(55, 85, 155, 85, 105, 180, 55, 85)

        canvas.pack(fill=BOTH, expand=1)

        container_frame = Frame(canvas, relief=RAISED, borderwidth=1)
        container_frame.pack(side='top', fill=BOTH, expand=True)

        for i in range(0,ncols):
            #primary.columnconfigure(i, pad=3)
            container_frame.columnconfigure(i, pad=3)

        for i in range(0,nrows+1):
            #primary.rowconfigure(i, pad=3)
            container_frame.rowconfigure(i, pad=3)

        #Fill a list with the actual buttons.
        #buttons = []
        for i in range(0,channels):
            self.button = Button(container_frame, text=str(pmt_mods[i])+'\n'+connections[str(i+1)][12], width=4, height=1)
            #button["command"] = self.fconnections
            self.button.bind("<Button-1>")#, self.fconnections)
            buttons.append(self.button)

        for i in range(0,nrows):
            for j in range(0,ncols):
                buttons[i*ncols+j].grid(row=i, column=j)

        #canvas.pack(fill=BOTH, expand=1)

def main():

    root = Tk()
    ex = Example()
    root.geometry("400x250+300+300")
    root.mainloop()


if __name__ == '__main__':
    main()
