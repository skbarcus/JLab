#!/usr/bin/env python3

"""
ZetCode Tkinter tutorial

In this script, we use the grid manager
to create a skeleton of a calculator.

Author: Jan Bodnar
Website: www.zetcode.com
"""

from tkinter import Tk, W, E
from tkinter.ttk import Frame, Button, Entry, Style

nrows = 24
ncols = 12
channels = nrows * ncols
pmt_mods = []    #Holds PMT module numbers.
buttons = []
#buttons = {}

#Dictionary containing all of the relevant connections for each PMT.
connections = {1:["b6-01","f6-00","L6.0"],2:["b6-03","f6-02","L7.0"],3:["b6-05","f6-04","L6.1"],4:["b6-07","f6-06","L7.1"],5:["b6-09","f6-08","L8.0"],6:["b6-11","f6-10","L9.0"],7:["b6-02","f6-01","L8.1"],8:["b6-04","f6-03","L9.1"],9:["b6-06","f6-05","L10.0"],10:["b6-08","f6-07","L11.0"],11:["b6-10","f6-09","L10.1"],12:["b6-12","f6-11","L11.1"],
               13:["b4-13","f4-12","L6.2"],14:["b4-15","f4-14","L7.2"],15:["b5-13","f5-12","L6.3."],16:["b5-15","f5-14","L7.3"],17:["b-","f-","L."],18:["b-","f-","L."],19:["b-","f-","L."],20:["b-","f-","L."],21:["b-","f-","L."],22:["b-","f-","L."],23:["b-","f-","L."],24:["b-","f-","L."],
               25:["b-","f-","L."],26:["b-","f-","L."],27:["b-","f-","L."],28:["b-","f-","L."],29:["b-","f-","L."],30:["b-","f-","L."],31:["b-","f-","L."],32:["b-","f-","L."],33:["b-","f-","L."],34:["b-","f-","L."],35:["b-","f-","L."],36:["b-","f-","L."],
               37:["b-","f-","L."],38:["b-","f-","L."],39:["b-","f-","L."],40:["b-","f-","L."],41:["b-","f-","L."],42:["b-","f-","L."],43:["b-","f-","L."],44:["b-","f-","L."],45:["b-","f-","L."],46:["b-","f-","L."],47:["b-","f-","L."],48:["b-","f-","L."],
               49:["b-","f-","L."],50:["b-","f-","L."],51:["b-","f-","L."],52:["b-","f-","L."],53:["b-","f-","L."],54:["b-","f-","L."],55:["b-","f-","L."],56:["b-","f-","L."],57:["b-","f-","L."],58:["b-","f-","L."],59:["b-","f-","L."],60:["b-","f-","L."],
               61:["b-","f-","L."],62:["b-","f-","L."],63:["b-","f-","L."],64:["b-","f-","L."],65:["b-","f-","L."],66:["b-","f-","L."],67:["b-","f-","L."],68:["b-","f-","L."],69:["b-","f-","L."],70:["b-","f-","L."],71:["b-","f-","L."],72:["b-","f-","L."],
               73:["b-","f-","L."],74:["b-","f-","L."],75:["b-","f-","L."],76:["b-","f-","L."],77:["b-","f-","L."],78:["b-","f-","L."],79:["b-","f-","L."],80:["b-","f-","L."],81:["b-","f-","L."],82:["b-","f-","L."],83:["b-","f-","L."],84:["b-","f-","L."],
               85:["b-","f-","L."],86:["b-","f-","L."],87:["b-","f-","L."],88:["b-","f-","L."],89:["b-","f-","L."],90:["b-","f-","L."],91:["b-","f-","L."],92:["b-","f-","L."],93:["b-","f-","L."],94:["b-","f-","L."],95:["b-","f-","L."],96:["b-","f-","L."],
               97:["b-","f-","L."],98:["b-","f-","L."],99:["b-","f-","L."],100:["b-","f-","L."],101:["b-","f-","L."],102:["b-","f-","L."],103:["b-","f-","L."],104:["b-","f-","L."],105:["b-","f-","L."],106:["b-","f-","L."],107:["b-","f-","L."],108:["b-","f-","L."],
               109:["b-","f-","L."],110:["b-","f-","L."],111:["b-","f-","L."],112:["b-","f-","L."],113:["b-","f-","L."],114:["b-","f-","L."],115:["b-","f-","L."],116:["b-","f-","L."],117:["b-","f-","L."],118:["b-","f-","L."],119:["b-","f-","L."],120:["b-","f-","L."],
               121:["b-","f-","L."],122:["b-","f-","L."],123:["b-","f-","L."],124:["b-","f-","L."],125:["b-","f-","L."],126:["b-","f-","L."],127:["b-","f-","L."],128:["b-","f-","L."],129:["b-","f-","L."],130:["b-","f-","L."],131:["b-","f-","L."],132:["b-","f-","L."],
               133:["b-","f-","L."],134:["b-","f-","L."],135:["b-","f-","L."],136:["b-","f-","L."],137:["b-","f-","L."],138:["b-","f-","L."],139:["b-","f-","L."],140:["b-","f-","L."],141:["b-","f-","L."],142:["b-","f-","L."],143:["b-","f-","L."],144:["b-","f-","L."],
               145:["a-","f-","L."],146:["a-","f-","L."],147:["a-","f-","L."],148:["a-","f-","L."],149:["a-","f-","L."],150:["a-","f-","L."],151:["a-","f-","L."],152:["a-","f-","L."],153:["a-","f-","L."],154:["a-","f-","L."],155:["a-","f-","L."],156:["a-","f-","L."],
               157:["a-","f-","L."],158:["a-","f-","L."],159:["a-","f-","L."],160:["a-","f-","L."],161:["a-","f-","L."],162:["a-","f-","L."],163:["a-","f-","L."],164:["a-","f-","L."],165:["a-","f-","L."],166:["a-","f-","L."],167:["a-","f-","L."],168:["a-","f-","L."],
               169:["a-","f-","L."],170:["a-","f-","L."],171:["a-","f-","L."],172:["a-","f-","L."],173:["a-","f-","L."],174:["a-","f-","L."],175:["a-","f-","L."],176:["a-","f-","L."],177:["a-","f-","L."],178:["a-","f-","L."],179:["a-","f-","L."],180:["a-","f-","L."],
               181:["a-","f-","L."],182:["a-","f-","L."],183:["a-","f-","L."],184:["a-","f-","L."],185:["a-","f-","L."],186:["a-","f-","L."],187:["a-","f-","L."],188:["a-","f-","L."],189:["a-","f-","L."],190:["a-","f-","L."],191:["a-","f-","L."],192:["a-","f-","L."],
               193:["a-","f-","L."],194:["a-","f-","L."],195:["a-","f-","L."],196:["a-","f-","L."],197:["a-","f-","L."],198:["a-","f-","L."],199:["a-","f-","L."],200:["a-","f-","L."],201:["a-","f-","L."],202:["a-","f-","L."],203:["a-","f-","L."],204:["a-","f-","L."],
               205:["a-","f-","L."],206:["a-","f-","L."],207:["a-","f-","L."],208:["a-","f-","L."],209:["a-","f-","L."],210:["a-","f-","L."],211:["a-","f-","L."],212:["a-","f-","L."],213:["a-","f-","L."],214:["a-","f-","L."],215:["a-","f-","L."],216:["a-","f-","L."],
               217:["a-","f-","L."],218:["a-","f-","L."],219:["a-","f-","L."],220:["a-","f-","L."],221:["a-","f-","L."],222:["a-","f-","L."],223:["a-","f-","L."],224:["a-","f-","L."],225:["a-","f-","L."],226:["a-","f-","L."],227:["a-","f-","L."],228:["a-","f-","L."],
               229:["a-","f-","L."],230:["a-","f-","L."],231:["a-","f-","L."],232:["a-","f-","L."],233:["a-","f-","L."],234:["a-","f-","L."],235:["a-","f-","L."],236:["a-","f-","L."],237:["a-","f-","L."],238:["a-","f-","L."],239:["a-","f-","L."],240:["a-","f-","L."],
               241:["a-","f-","L."],242:["a-","f-","L."],243:["a-","f-","L."],244:["a-","f-","L."],245:["a-","f-","L."],246:["a-","f-","L."],247:["a-","f-","L."],248:["a-","f-","L."],249:["a-","f-","L."],250:["a-","f-","L."],251:["a-","f-","L."],252:["a-","f-","L."],
               253:["a-","f-","L."],254:["a-","f-","L."],255:["a-","f-","L."],256:["a-","f-","L."],257:["a-","f-","L."],258:["a-","f-","L."],259:["a-","f-","L."],260:["a-","f-","L."],261:["a-","f-","L."],262:["a-","f-","L."],263:["a-","f-","L."],264:["a-","f-","L."],
               265:["a-","f-","L."],266:["a-","f-","L."],267:["a-","f-","L."],268:["a-","f-","L."],269:["a-","f-","L."],270:["a-","f-","L."],271:["a-","f-","L."],272:["a-","f-","L."],273:["a-","f-","L."],274:["a-","f-","L."],275:["a-","f-","L."],276:["a-","f-","L."],
               277:["a-","f-","L."],278:["a-","f-","L."],279:["a-","f-","L."],280:["a-","f-","L."],281:["a-","f-","L."],282:["a-","f-","L."],283:["a-","f-","L."],284:["a-","f-","L."],285:["a-","f-","L."],286:["a-","f-","L."],287:["a-","f-","L."],288:["a-","f-","L."]}
#print(connections)

#Fill PMT module numbers.
for i in range(0,channels):
    pmt_mods.append(i+1)
    #print(pmt_mods[i])

#Fill a list with the names for the buttons.
button_names = []
for i in range(0,channels):
    button_names.append("Button"+str(i+1))
#print(button_names)

class Test(Frame):

    def __init__(self):
        super().__init__()
        self.initUI()
        global buttons
        #self.fconnections()


    def initUI(self):
        global buttons
        self.master.title("HCal Detector")

        Style().configure("TButton", padding=(0, 5, 0, 5),
            font='serif 10')

        #entry = Entry(self)
        #entry.grid(row=0, columnspan=4, sticky=W+E)

        for i in range(0,ncols):
            self.columnconfigure(i, pad=3)

        for i in range(0,nrows):
            self.rowconfigure(i, pad=3)

        #Fill a list with the actual buttons.
        #buttons = []
        for i in range(0,channels):
            button = Button(self, text=pmt_mods[i])#, command=self.fconnections)
            #button = Button(self, text=pmt_mods[i], command=lambda event, obj=button: fconnections(event,obj))
            button["command"] = self.fconnections
            buttons.append(button)
            #buttons.append(Button(self, text=pmt_mods[i], command=self.fconnections))
            #buttons[i].bind("<Button-1>", self.fconnections)
            #button = Button(self, text=pmt_mods[i], command=self.fconnections)
            #buttons[button] = pmt_mods[i]
            #print(button_names[i]," = ",buttons[i])

        for i in range(0,nrows):
            for j in range(0,ncols):
                buttons[i*ncols+j].grid(row=i, column=j)

        self.pack()

    def fconnections(self,event):
        global buttons
        print('PMT Module ',i,': amplifier channel = ', connections[i][0],', fADC channel = ', connections[i][1], ', HV channel = ',connections[i][2],'.')
        #print(buttons[0].cget('text'))
        #print(event.widget.cget('text'))
        #print(event.cget('text'))
        #print(example.event.widget)
        print(buttons[10])
        print(event.widget)
        #print("Hello there!")

    #def connections(event,obj):
        #print('clicked ',obj)

def main():

    root = Tk()
    app = Test()
    root.mainloop()


if __name__ == '__main__':
    main()
