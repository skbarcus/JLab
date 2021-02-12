from tkinter import Tk, Label, Button, StringVar, W, E, RIGHT, BOTH, RAISED, LEFT, CENTER, TOP, BOTTOM, Canvas, Scrollbar
#import tkinter.ttk as ttk
from tkinter.ttk import Style, Entry, Frame
#import tkinter.font as tkFont

nrows = 24
ncols = 12
channels = nrows * ncols
pmt_mods = []    #Holds PMT module numbers.
buttons = []

#Dictionary containing all of the relevant connections for each PMT.
connections = {1:["b6-01","f6-00","L6.0"],2:["b6-03","f6-02","L7.0"],3:["b6-05","f6-04","L6.1"],4:["b6-07","f6-06","L7.1"],5:["b6-09","f6-08","L8.0"],6:["b6-11","f6-10","L9.0"],7:["b6-02","f6-01","L8.1"],8:["b6-04","f6-03","L9.1"],9:["b6-06","f6-05","L10.0"],10:["b6-08","f6-07","L11.0"],11:["b6-10","f6-09","L10.1"],12:["b6-12","f6-11","L11.1"],
               13:["b4-13","f4-12","L6.2"],14:["b4-15","f4-14","L7.2"],15:["b5-13","f5-12","L6.3"],16:["b5-15","f5-14","L7.3"],17:["b6-13","f6-12","L8.2"],18:["b6-15","f6-14","L9.2"],19:["b4-14","f4-13","L8.3"],20:["b4-16","f4-15","L9.3"],21:["b5-14","f5-13","L10.2"],22:["b5-16","f5-15","L11.2"],23:["b6-14","f6-13","L10.3"],24:["b6-16","f6-15","L11.3"],
               25:["b7-01","f7-00","L6.4"],26:["b7-03","f7-02","L7.4"],27:["b7-05","f7-04","L6.5"],28:["b7-07","f7-06","L7.5"],29:["b7-09","f7-08","L8.4"],30:["b7-11","f7-10","L9.4"],31:["b7-02","f7-01","L8.5"],32:["b7-04","f7-03","L9.5"],33:["b7-06","f7-05","L10.4"],34:["b7-08","f7-07","L11.4"],35:["b7-10","f7-09","L10.5"],36:["b7-12","f7-11","L11.5"],
               37:["b8-01","f8-00","L6.6"],38:["b8-03","f8-02","L7.6"],39:["b8-05","f8-04","L6.7"],40:["b8-07","f8-06","L7.7"],41:["b8-09","f8-08","L8.6"],42:["b8-11","f8-10","L9.6"],43:["b8-02","f8-01","L8.7"],44:["b8-04","f8-03","L9.7"],45:["b8-06","f8-05","L10.6"],46:["b8-08","f8-07","L11.6"],47:["b8-10","f8-09","L10.7"],48:["b8-12","f8-11","L11.7"],
               49:["b9-01","f9-00","L6.8"],50:["b9-03","f9-02","L7.8"],51:["b9-05","f9-04","L6.9"],52:["b9-07","f9-06","L7.9"],53:["b9-09","f9-08","L8.8"],54:["b9-11","f9-10","L9.8"],55:["b9-02","f9-01","L8.9"],56:["b9-04","f9-03","L9.9"],57:["b9-06","f9-05","L10.8"],58:["b9-08","f9-07","L11.8"],59:["b9-10","f9-09","L10.9"],60:["b9-12","f9-11","L11.9"],
               61:["b7-13","f7-12","L6.10"],62:["b7-15","f7-14","L7.10"],63:["b8-13","f8-12","L6.11"],64:["b8-15","f8-14","L7.11"],65:["b9-13","f9-12","L8.10"],66:["b9-15","f9-14","L9.10"],67:["b7-14","f7-13","L8.11"],68:["b7-16","f7-15","L9.11"],69:["b8-14","f8-13","L10.10"],70:["b8-16","f8-15","L11.10"],71:["b9-14","f9-13","L10.11"],72:["b9-16","f9-15","L11.11"],
               73:["b1-01","f1-00","L0.0"],74:["b1-03","f1-02","L1.0"],75:["b1-05","f1-04","L0.1"],76:["b1-07","f1-06","L1.1"],77:["b1-09","f1-08","L2.0"],78:["b1-11","f1-10","L3.0"],79:["b1-02","f1-01","L2.1"],80:["b1-04","f1-03","L3.1"],81:["b1-06","f1-05","L4.0"],82:["b1-08","f1-07","L5.0"],83:["b1-10","f1-09","L4.1"],84:["b1-12","f1-11","L5.1"],
               85:["b2-01","f2-00","L0.2"],86:["b2-03","f2-02","L1.2"],87:["b2-05","f2-04","L0.3"],88:["b2-09","f2-08","L1.3"],89:["b2-09","f2-08","L2.2"],90:["b2-11","f2-10","L3.2"],91:["b2-02","f2-01","L2.3"],92:["b2-04","f2-03","L3.3"],93:["b2-06","f2-05","L4.2"],94:["b2-08","f2-07","L5.2"],95:["b2-10","f2-09","L4.3"],96:["b2-12","f2-11","L5.3"],
               97:["b3-01","f3-00","L0.4"],98:["b3-03","f3-02","L1.4"],99:["b3-05","f3-04","L0.5"],100:["b3-07","f3-06","L1.5"],101:["b3-09","f3-08","L2.4"],102:["b3-11","f3-10","L3.4"],103:["b3-02","f3-01","L2.5"],104:["b3-04","f3-03","L3.5"],105:["b3-06","f3-05","L4.4"],106:["b3-08","f3-07","L5.4"],107:["b3-10","f3-09","L4.5"],108:["b3-12","f3-11","L5.5"],
               109:["b1-13","f1-12","L0.6"],110:["b1-15","f1-14","L1.6"],111:["b2-13","f2-12","L0.7"],112:["b2-15","f2-14","L1.7"],113:["b3-13","f3-12","L2.6"],114:["b3-15","f3-14","L3.6"],115:["b1-14","f1-13","L2.7"],116:["b1-16","f1-15","L3.7"],117:["b2-14","f2-13","L4.6"],118:["b2-16","f2-15","L5.6"],119:["b3-14","f3-13","L4.7"],120:["b3-16","f3-15","L5.7"],
               121:["b4-01","f4-00","L0.8"],122:["b4-03","f4-02","L1.8"],123:["b4-05","f4-04","L0.9"],124:["b4-07","f4-06","L1.9"],125:["b4-09","f4-08","L2.8"],126:["b4-11","f4-10","L3.8"],127:["b4-02","f4-01","L2.9"],128:["b4-04","f4-03","L3.9"],129:["b4-06","f4-05","L4.8"],130:["b4-08","f4-07","L5.8"],131:["b4-10","f4-09","L4.9"],132:["b4-12","f4-11","L5.9"],
               133:["b5-01","f5-00","L0.10"],134:["b5-03","f5-02","L1.10"],135:["b5-05","f5-04","L0.11"],136:["b5-07","f5-06","L1.11"],137:["b5-09","f5-08","L2.10"],138:["b5-11","f5-10","L3.10"],139:["b5-02","f5-01","L2.11"],140:["b5-04","f5-03","L3.11"],141:["b5-06","f5-05","L4.10"],142:["b5-08","f5-07","L5.10"],143:["b5-10","f5-09","L4.11"],144:["b5-12","f5-11","L5.11"],
               145:["a1-01","f10-00","L0.0"],146:["a1-03","f10-02","L1.0"],147:["a1-05","f10-04","L0.1"],148:["a1-07","f10-06","L1.1"],149:["a1-09","f10-08","L3.0"],150:["a1-11","f10-10","L4.0"],151:["a1-02","f10-01","L3.1"],152:["a1-04","f10-03","L4.1"],153:["a1-06","f10-05","L5.0"],154:["a1-08","f10-07","L6.0"],155:["a1-10","f10-09","L5.1"],156:["a1-12","f10-11","L6.1"],
               157:["a2-01","f11-00","L0.2"],158:["a2-03","f11-02","L1.2"],159:["a2-05","f11-04","L0.3"],160:["a2-07","f11-06","L1.3"],161:["a2-09","f11-08","L3.2"],162:["a2-11","f11-10","L4.2"],163:["a2-02","f11-01","L3.3"],164:["a2-04","f11-03","L4.3"],165:["a2-06","f11-05","L5.2"],166:["a2-08","f11-07","L6.2"],167:["a2-10","f11-09","L5.3"],168:["a2-12","f11-11","L6.3"],
               169:["a3-01","f12-00","L0.4"],170:["a3-03","f12-02","L1.4"],171:["a3-05","f12-04","L0.5"],172:["a3-07","f12-06","L1.5"],173:["a3-09","f12-08","L3.4"],174:["a3-11","f12-10","L4.4"],175:["a3-02","f12-01","L3.5"],176:["a3-04","f12-03","L4.5"],177:["a3-06","f12-05","L5.4"],178:["a3-08","f12-07","L6.4"],179:["a3-10","f12-09","L5.5"],180:["a3-12","f12-11","L6.5"],
               181:["a1-13","f10-12","L0.6"],182:["a1-15","f10-14","L1.6"],183:["a2-13","f11-12","L0.7"],184:["a2-15","f11-14","L1.7"],185:["a3-13","f12-12","L3.6"],186:["a3-15","f12-14","L4.6"],187:["a1-14","f10-13","L3.7"],188:["a1-16","f10-15","L4.7"],189:["a2-14","f11-13","L5.6"],190:["a2-16","f11-15","L6.6"],191:["a3-14","f12-13","L5.7"],192:["a3-16","f12-15","L6.7"],
               193:["a4-01","f13-00","L0.8"],194:["a4-03","f13-02","L1.8"],195:["a4-05","f13-04","L0.9"],196:["a4-07","f13-06","L1.9"],197:["a4-09","f13-08","L3.8"],198:["a4-11","f13-10","L4.8"],199:["a4-02","f13-01","L3.9"],200:["a4-04","f13-03","L4.9"],201:["a4-06","f13-05","L5.8"],202:["a4-08","f13-07","L6.8"],203:["a4-10","f13-09","L5.9"],204:["a4-12","f13-11","L6.9"],
               205:["a5-01","f14-00","L0.10"],206:["a5-03","f14-02","L1.10"],207:["a5-05","f14-04","L0.11"],208:["a5-07","f14-06","L1.11"],209:["a5-09","f14-08","L3.10"],210:["a5-11","f14-10","L4.10"],211:["a5-02","f14-01","L3.11"],212:["a5-04","f14-03","L4.11"],213:["a5-06","f14-05","L5.10"],214:["a5-08","f14-07","L6.10"],215:["a5-10","f14-09","L5.11"],216:["a5-12","f14-11","L6.11"],
               217:["a6-01","f15-00","L7.0"],218:["a6-03","f15-02","L8.0"],219:["a6-05","f15-04","L7.1"],220:["a6-07","f15-06","L8.1"],221:["a6-09","f15-08","L9.0"],222:["a6-11","f15-10","L10.0"],223:["a6-02","f15-01","L9.1"],224:["a6-04","f15-03","L10.1"],225:["a6-06","f15-05","L11.0"],226:["a6-08","f15-07","L12.0."],227:["a6-10","f15-09","L11.1"],228:["a6-12","f15-11","L12.1"],
               229:["a4-13","f13-12","L7.2"],230:["a4-15","f13-14","L8.2"],231:["a5-13","f14-12","L7.3"],232:["a5-15","f14-14","L8.3"],233:["a6-13","f15-12","L9.2"],234:["a6-15","f15-14","L10.2"],235:["a4-14","f13-13","L9.3"],236:["a4-16","f13-15","L10.3"],237:["a5-14","f14-13","L11.2"],238:["a5-16","f14-15","L12.2"],239:["a6-14","f15-13","L11.3"],240:["a6-16","f15-15","L12.3"],
               241:["a7-01","f16-00","L7.4"],242:["a7-03","f16-02","L8.4"],243:["a7-05","f16-04","L7.5"],244:["a7-07","f16-06","L8.5"],245:["a7-09","f16-08","L9.4"],246:["a7-11","f16-10","L10.4"],247:["a7-02","f16-01","L9.5"],248:["a7-04","f16-03","L10.5"],249:["a7-06","f16-05","L11.4"],250:["a7-08","f16-07","L12.4"],251:["a7-10","f16-09","L11.5"],252:["a7-12","f16-11","L12.5"],
               253:["a8-01","f17-00","L7.6"],254:["a8-03","f17-02","L8.6"],255:["a8-05","f17-04","L7.7"],256:["a8-07","f17-06","L8.7"],257:["a8-09","f17-08","L9.6"],258:["a8-11","f17-10","L10.6"],259:["a8-02","f17-01","L9.7"],260:["a8-04","f17-03","L10.7"],261:["a8-06","f17-05","L11.6"],262:["a8-08","f17-07","L12.6"],263:["a8-10","f17-09","L11.7"],264:["a8-12","f17-11","L12.7"],
               265:["a9-01","f18-00","L7.8"],266:["a9-03","f18-02","L8.8"],267:["a9-05","f18-04","L7.9"],268:["a9-07","f18-06","L8.9"],269:["a9-09","f18-08","L9.8"],270:["a9-11","f18-10","L10.8"],271:["a9-02","f18-01","L9.9"],272:["a9-04","f18-03","L10.9"],273:["a9-06","f18-05","L11.8"],274:["a9-08","f18-07","L12.8"],275:["a9-10","f18-09","L11.9"],276:["a9-12","f18-11","L12.9"],
               277:["a7-13","f16-12","L7.10"],278:["a7-15","f16-14","L8.10"],279:["a8-13","f17-12","L7.11"],280:["a8-15","f17-14","L8.11"],281:["a9-13","f18-12","L9.10"],282:["a9-15","f18-14","L10.10"],283:["a7-14","f16-13","L9.11"],284:["a7-16","f16-15","L10.11"],285:["a8-14","f17-13","L11.10"],286:["a8-16","f17-15","L12.10"],287:["a9-14","f18-13","L11.11"],288:["a9-16","f18-15","L12.11"]}
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

class MyFirstGUI:
    def __init__(self, primary):
        self.primary = primary
        primary.title("Front-End Electronics")

        self.style = Style()
        self.style.theme_use("default")

        #Create three frames to hold the contents of each FE rack.
        RR1_frame = Frame(primary, relief=RAISED, borderwidth=2)
        RR1_frame.pack(side='left', fill=BOTH, expand=True)

        RR2_frame = Frame(primary, relief=RAISED, borderwidth=2)
        RR2_frame.pack(side='left', fill=BOTH, expand=True)

        RR3_frame = Frame(primary, relief=RAISED, borderwidth=2)
        RR3_frame.pack(side='left', fill=BOTH, expand=True)

        #Create buttons and labels for electronics in RR1
        RR1_label = Label(RR1_frame, text='RR1', font='Helvetica 16 bold')
        RR1_label.pack(side="top")

        RR1_sum_mods_btn = Button(RR1_frame, text='Summing Modules', width=10, height=5, font='Helvetica 14 bold')
        RR1_sum_mods_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR1_sum_mods_btn.bind("<Button-1>", self.RR1_sum_connections)

        RR1_splitter_btn = Button(RR1_frame, text='Splitter Panels', width=10, height=13, font='Helvetica 14 bold')
        RR1_splitter_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR1_splitter_btn.bind("<Button-1>", self.RR1_splitter_connections)

        RR1_amp_btn = Button(RR1_frame, text='Amplifiers', width=10, height=5, font='Helvetica 14 bold')
        RR1_amp_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR1_amp_btn.bind("<Button-1>", self.RR1_amp_connections)

        #Create buttons and labels for electronics in RR2
        RR2_label = Label(RR2_frame, text='RR2', font='Helvetica 16 bold')
        RR2_label.pack(side="top")

        RR2_upper_tdc_pp_btn = Button(RR2_frame, text='Upper TDC Patch Panels', width=10, height=4, font='Helvetica 14 bold')
        RR2_upper_tdc_pp_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR2_upper_tdc_pp_btn.bind("<Button-1>", self.RR2_upper_tdc_pp_connections)

        RR2_upper_tdc_disc_btn = Button(RR2_frame, text='TDC Discriminators from RR3', width=10, height=5, font='Helvetica 14 bold')
        RR2_upper_tdc_disc_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR2_upper_tdc_disc_btn.bind("<Button-1>", self.RR2_Upper_TDC_Disc_connections)

        RR2_lower_tdc_pp_btn = Button(RR2_frame, text='Lower TDC Patch Panels', width=10, height=6, font='Helvetica 14 bold')
        RR2_lower_tdc_pp_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR2_lower_tdc_pp_btn.bind("<Button-1>", self.RR2_lower_tdc_pp_connections)

        RR2_lower_tdc_disc_btn = Button(RR2_frame, text='TDC Discriminators from RR1', width=10, height=5, font='Helvetica 14 bold')
        RR2_lower_tdc_disc_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR2_lower_tdc_disc_btn.bind("<Button-1>", self.RR2_Lower_TDC_Disc_connections)

        RR2_fadc_pp_btn = Button(RR2_frame, text='fADC Patch Panels', width=10, height=11, font='Helvetica 14 bold')
        RR2_fadc_pp_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR2_fadc_pp_btn.bind("<Button-1>", self.RR2_adc_pp_connections)

        #Create buttons and labels for electronics in RR3
        RR3_label = Label(RR3_frame, text='RR3', font='Helvetica 16 bold')
        RR3_label.pack(side="top")

        RR3_sum_mods_btn = Button(RR3_frame, text='Summing Modules', width=10, height=5, font='Helvetica 14 bold')
        RR3_sum_mods_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR3_sum_mods_btn.bind("<Button-1>", self.RR3_sum_connections)

        RR3_splitter_btn = Button(RR3_frame, text='Splitter Panels', width=10, height=13, font='Helvetica 14 bold')
        RR3_splitter_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR3_splitter_btn.bind("<Button-1>", self.RR3_splitter_connections)

        RR3_amp_btn = Button(RR3_frame, text='Amplifiers', width=10, height=5, font='Helvetica 14 bold')
        RR3_amp_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR3_amp_btn.bind("<Button-1>", self.RR3_amp_connections)

    def fconnections(self, event):
        pmt = event.widget.cget('text')
        print('PMT Module ',pmt,': amplifier channel = ', connections[pmt][0],', fADC channel = ', connections[pmt][1], ', HV channel = ',connections[pmt][2])
        #print(event.widget.cget('text'))

    def RR3_amp_connections(self, event):
        win = Tk()
        win.wm_title("RR3 Amplifier Connections")
        win.geometry("1400x600")

        RR3_amp_frames = []                           #List to hold the amplifier frames.
        buttons = []                                  #List to hold the amplifier buttons.
        nslots = 9                                    #Number of amplifiers in crate.
        nrows = 16                                    #Number of rows of channels in an amp.
        ncols = 3                                     #Number of columns of channels in an amp.
        nchannels = nrows * ncols                     #Number of channels in one amp.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the amps.

        #Create 9 different frames to hold the amp channel buttons.
        for i in range (0,nslots):
            frame = Frame(win, relief=RAISED, borderwidth=2)
            frame.pack(side='left', fill=BOTH, expand=True)
            RR3_amp_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            for j in range(0,ncols):
                RR3_amp_frames[i].columnconfigure(j, pad=3, weight=1)

            for j in range(0,nrows+1):
                RR3_amp_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the amp buttons to their functions.
        for i in range(0,nslots):
            buttons.append([])#Make this list 2D to hold the buttons of each of the 9 frames separately.
            for j in range(0,nrows):
                for k in range(0,ncols):

                    btn = Button(RR3_amp_frames[i], width=1, height=1, font='Helvetica 8')

                    #Give names and function binds to the three different columns.
                    if k == 0:
                        text = 'a'+str(i+1)+'-'+str(j+1)+'\n In'
                        btn['text'] = text
                        #btn['font'] = 'Helvetica 20'
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)
                    elif k == 1:
                        text = 'a'+str(i+1)+'-'+str(j+1)+'\n Out 1'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)
                    elif k == 2:
                        text = 'a'+str(i+1)+'-'+str(j+1)+'\n Out 2'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)

        labels = []
        #Place labels above the amp buttons in the grid.
        for i in range(0,nslots):
            label = Label(RR3_amp_frames[i], text='Amplifier '+str(i+10))
            labels.append(label)
            labels[i].grid(row=0,columnspan=ncols)

        #Place the amp buttons in the grid.
        for i in range(0,nslots):
            for j in range(1,nrows+1):
                for k in range(0,ncols):
                    buttons[i][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')

    def ftest(self, event):
        print('Test successful!', event.widget.cget('text'))

    def RR1_amp_connections(self, event):
        win = Tk()
        win.wm_title("RR1 Amplifier Connections")
        win.geometry("1400x600")

        RR1_amp_frames = []                           #List to hold the amplifier frames.
        buttons = []                                  #List to hold the amplifier buttons.
        nslots = 9                                    #Number of amplifiers in crate.
        nrows = 16                                    #Number of rows of channels in an amp.
        ncols = 3                                     #Number of columns of channels in an amp.
        nchannels = nrows * ncols                     #Number of channels in one amp.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the amps.

        #Create 9 different frames to hold the amp channel buttons.
        for i in range (0,nslots):
            frame = Frame(win, relief=RAISED, borderwidth=2)
            frame.pack(side='left', fill=BOTH, expand=True)
            RR1_amp_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            for j in range(0,ncols):
                RR1_amp_frames[i].columnconfigure(j, pad=3, weight=1)

            for j in range(0,nrows+1):
                RR1_amp_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the amp buttons to their functions.
        for i in range(0,nslots):
            buttons.append([])#Make this list 2D to hold the buttons of each of the 9 frames separately.
            for j in range(0,nrows):
                for k in range(0,ncols):

                    btn = Button(RR1_amp_frames[i], width=1, height=1, font='Helvetica 8')

                    #Give names and function binds to the three different columns.
                    if k == 0:
                        text = 'b'+str(i+1)+'-'+str(j+1)+'\n In'
                        btn['text'] = text
                        #btn['font'] = 'Helvetica 20'
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)
                    elif k == 1:
                        text = 'b'+str(i+1)+'-'+str(j+1)+'\n Out 1'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)
                    elif k == 2:
                        text = 'b'+str(i+1)+'-'+str(j+1)+'\n Out 2'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)

        labels = []
        #Place labels above the amp buttons in the grid.
        for i in range(0,nslots):
            label = Label(RR1_amp_frames[i], text='Amplifier '+str(i+1))
            labels.append(label)
            labels[i].grid(row=0,columnspan=ncols)

        #Place the amp buttons in the grid.
        for i in range(0,nslots):
            for j in range(1,nrows+1):
                for k in range(0,ncols):
                    buttons[i][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')

    def RR3_sum_connections(self, event):
        win = Tk()
        win.wm_title("RR3 Summing Module Connections")
        win.geometry("1400x600")

        RR3_sum_frames = []                           #List to hold the frames.
        buttons = []                                  #List to hold the  buttons.
        nslots = 12                                   #Number of slots in crate.
        nsum_mods = 5                                 #Number of summing modules.
        nrows = 16                                    #Number of rows of channels in a module.
        ncols = 2                                     #Number of columns of channels in a module.
        nchannels = nrows * ncols                     #Number of channels in one module.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the modules.

        #Create 9 different frames to hold the sum channel buttons.
        for i in range (0,nslots):
            frame = Frame(win, relief=RAISED, borderwidth=2)
            frame.pack(side='left', fill=BOTH, expand=True)
            RR3_sum_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            if i == 0 or i == 2 or i == 4 or i == 6 or i == 8:
                for j in range(0,ncols):
                    RR3_sum_frames[i].columnconfigure(j, pad=3, weight=1)

                for j in range(0,nrows+5):
                    RR3_sum_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the sum buttons to their functions.
        nmods = 0
        for i in range(0,nslots):
            if i == 0 or i == 2 or i == 4 or i == 6 or i == 8:
                buttons.append([])#Make this list 2D to hold the buttons of each of the frames separately.
                for j in range(0,nrows+4):
                    for k in range(0,ncols):
                        btn = Button(RR3_sum_frames[i], width=1, height=1, font='Helvetica 8')

                        #Give names and function binds to the three different columns.
                        if j<nrows:
                            if k == 0:
                                text = str(nmods+6)+'A-'+str(j+1)+'\n In'
                                btn['text'] = text
                                btn.bind("<Button-1>", self.ftest)
                                buttons[nmods].append(btn)
                            elif k == 1:
                                text = str(nmods+6)+'B-'+str(j+1)+'\n In'
                                btn['text'] = text
                                btn.bind("<Button-1>", self.ftest)
                                buttons[nmods].append(btn)
                        #Create different buttons for the outputs.
                        else:
                            if k == 0:
                                text = str(nmods+6)+'A \n Out'+str(j-nrows)
                                btn['text'] = text
                                btn.bind("<Button-1>", self.ftest)
                                buttons[nmods].append(btn)
                            elif k == 1:
                                text = str(nmods+6)+'B \n Out'+str(j-nrows)
                                btn['text'] = text
                                btn.bind("<Button-1>", self.ftest)
                                buttons[nmods].append(btn)
                nmods = nmods + 1

        #Place labels above the sum buttons in the grid.
        labels = []
        nlabels = 0
        for i in range(0,nslots):
            if i == 0 or i == 2 or i == 4 or i == 6 or i == 8:
                label = Label(RR3_sum_frames[i], text='Summing Module '+str(nlabels+6))
                labels.append(label)
                labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                    
        #Place the sum buttons in the grid.
        nmods = 0
        for i in range(0,nslots):
            if i == 0 or i == 2 or i == 4 or i == 6 or i == 8:
                for j in range(1,nrows+5):
                    for k in range(0,ncols):
                        buttons[nmods][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')
                nmods = nmods +1

    def RR1_sum_connections(self, event):
        win = Tk()
        win.wm_title("RR1 Summing Module Connections")
        win.geometry("1400x600")

        RR1_sum_frames = []                           #List to hold the frames.
        buttons = []                                  #List to hold the  buttons.
        nslots = 12                                   #Number of slots in crate.
        nsum_mods = 5                                 #Number of summing modules.
        nrows = 16                                    #Number of rows of channels in a module.
        ncols = 2                                     #Number of columns of channels in a module.
        nchannels = nrows * ncols                     #Number of channels in one module.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the modules.

        #Create 9 different frames to hold the sum channel buttons.
        for i in range (0,nslots):
            frame = Frame(win, relief=RAISED, borderwidth=2)
            frame.pack(side='left', fill=BOTH, expand=True)
            RR1_sum_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            if i == 0 or i == 2 or i == 4 or i == 6 or i == 8:
                for j in range(0,ncols):
                    RR1_sum_frames[i].columnconfigure(j, pad=3, weight=1)

                for j in range(0,nrows+5):
                    RR1_sum_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the sum buttons to their functions.
        nmods = 0
        for i in range(0,nslots):
            if i == 0 or i == 2 or i == 4 or i == 6 or i == 8:
                buttons.append([])#Make this list 2D to hold the buttons of each of the frames separately.
                for j in range(0,nrows+4):
                    for k in range(0,ncols):
                        btn = Button(RR1_sum_frames[i], width=1, height=1, font='Helvetica 8')

                        #Give names and function binds to the three different columns.
                        if j<nrows:
                            if k == 0:
                                text = str(nmods+1)+'A-'+str(j+1)+'\n In'
                                btn['text'] = text
                                btn.bind("<Button-1>", self.ftest)
                                buttons[nmods].append(btn)
                            elif k == 1:
                                text = str(nmods+1)+'B-'+str(j+1)+'\n In'
                                btn['text'] = text
                                btn.bind("<Button-1>", self.ftest)
                                buttons[nmods].append(btn)
                        #Create different buttons for the outputs.
                        else:
                            if k == 0:
                                text = str(nmods+1)+'A \n Out'+str(j-nrows)
                                btn['text'] = text
                                btn.bind("<Button-1>", self.ftest)
                                buttons[nmods].append(btn)
                            elif k == 1:
                                text = str(nmods+1)+'B \n Out'+str(j-nrows)
                                btn['text'] = text
                                btn.bind("<Button-1>", self.ftest)
                                buttons[nmods].append(btn)
                nmods = nmods + 1

        #Place labels above the sum buttons in the grid.
        labels = []
        nlabels = 0
        for i in range(0,nslots):
            if i == 0 or i == 2 or i == 4 or i == 6 or i == 8:
                label = Label(RR1_sum_frames[i], text='Summing Module '+str(nlabels+1))
                labels.append(label)
                labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                    
        #Place the sum buttons in the grid.
        nmods = 0
        for i in range(0,nslots):
            if i == 0 or i == 2 or i == 4 or i == 6 or i == 8:
                for j in range(1,nrows+5):
                    for k in range(0,ncols):
                        buttons[nmods][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')
                nmods = nmods +1

    def RR2_Upper_TDC_Disc_connections(self, event):
        win = Tk()
        win.wm_title("RR2 Upper TDC Discriminator Connections")
        win.geometry("1400x600")
        RR2_upper_TDC_disc_frames = []                #List to hold the frames.
        buttons = []                                  #List to hold the  buttons.
        nslots = 12                                   #Number of slots in crate.
        ndisc = 9                                     #Number of summing modules.
        nrows = 16                                    #Number of rows of channels in a module.
        ncols = 3                                     #Number of columns of channels in a module.
        nchannels = nrows * ncols                     #Number of channels in one module.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the modules.

        #Create different frames to hold the channel buttons.
        for i in range (0,nslots):
            frame = Frame(win, relief=RAISED, borderwidth=2)
            frame.pack(side='left', fill=BOTH, expand=True)
            RR2_upper_TDC_disc_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            if i == 0 or i == 1 or i == 3 or i == 4 or i == 6 or i == 7 or i == 8 or i == 10 or i == 11:
                for j in range(0,ncols):
                    RR2_upper_TDC_disc_frames[i].columnconfigure(j, pad=3, weight=1)

                for j in range(0,nrows+1):
                    RR2_upper_TDC_disc_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the sum buttons to their functions.
        nmods = 0
        for i in range(0,nslots):
            if i == 0 or i == 1 or i == 3 or i == 4 or i == 6 or i == 7 or i == 8 or i == 10 or i == 11:
                buttons.append([])#Make this list 2D to hold the buttons of each of the frames separately.
                for j in range(0,nrows):
                    for k in range(0,ncols):
                        btn = Button(RR2_upper_TDC_disc_frames[i], width=1, height=1, font='Helvetica 8')

                        #Give names and function binds to the three different columns.
                        if k == 0:
                            text = str(18-nmods)+'-'+str(j+1)+'\n In'
                            btn['text'] = text
                            btn.bind("<Button-1>", self.ftest)
                            buttons[nmods].append(btn)
                        elif k == 1:
                            text = str(18-nmods)+'-'+str(j+1)+'\n Out 1'
                            btn['text'] = text
                            btn.bind("<Button-1>", self.ftest)
                            buttons[nmods].append(btn)
                        elif k == 2:
                            text = str(18-nmods)+'-'+str(j+1)+'\n Out 2'
                            btn['text'] = text
                            btn.bind("<Button-1>", self.ftest)
                            buttons[nmods].append(btn)
                nmods = nmods + 1

        #Place labels above the sum buttons in the grid.
        labels = []
        nlabels = 0
        for i in range(0,nslots):
            if i == 0 or i == 1 or i == 3 or i == 4 or i == 6 or i == 7 or i == 8 or i == 10 or i == 11:
                label = Label(RR2_upper_TDC_disc_frames[i], text='TDC Discriminator '+str(18-nlabels))
                labels.append(label)
                labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                    
        #Place the sum buttons in the grid.
        nmods = 0
        for i in range(0,nslots):
            if i == 0 or i == 1 or i == 3 or i == 4 or i == 6 or i == 7 or i == 8 or i == 10 or i == 11:
                for j in range(1,nrows+1):
                    for k in range(0,ncols):
                        buttons[nmods][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')
                nmods = nmods +1

    def RR2_Lower_TDC_Disc_connections(self, event):
        win = Tk()
        win.wm_title("RR2 Lower TDC Discriminator Connections")
        win.geometry("1400x600")
        RR2_lower_TDC_disc_frames = []                #List to hold the frames.
        buttons = []                                  #List to hold the  buttons.
        nslots = 12                                   #Number of slots in crate.
        ndisc = 9                                     #Number of summing modules.
        nrows = 16                                    #Number of rows of channels in a module.
        ncols = 3                                     #Number of columns of channels in a module.
        nchannels = nrows * ncols                     #Number of channels in one module.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the modules.

        #Create different frames to hold the channel buttons.
        for i in range (0,nslots):
            frame = Frame(win, relief=RAISED, borderwidth=2)
            frame.pack(side='left', fill=BOTH, expand=True)
            RR2_lower_TDC_disc_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            if i == 0 or i == 1 or i == 3 or i == 4 or i == 6 or i == 7 or i == 8 or i == 10 or i == 11:
                for j in range(0,ncols):
                    RR2_lower_TDC_disc_frames[i].columnconfigure(j, pad=3, weight=1)

                for j in range(0,nrows+1):
                    RR2_lower_TDC_disc_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the sum buttons to their functions.
        nmods = 0
        for i in range(0,nslots):
            if i == 0 or i == 1 or i == 3 or i == 4 or i == 6 or i == 7 or i == 8 or i == 10 or i == 11:
                buttons.append([])#Make this list 2D to hold the buttons of each of the frames separately.
                for j in range(0,nrows):
                    for k in range(0,ncols):
                        btn = Button(RR2_lower_TDC_disc_frames[i], width=1, height=1, font='Helvetica 8')

                        #Give names and function binds to the three different columns.
                        if k == 0:
                            text = str(nmods+1)+'-'+str(j+1)+'\n In'
                            btn['text'] = text
                            btn.bind("<Button-1>", self.ftest)
                            buttons[nmods].append(btn)
                        elif k == 1:
                            text = str(nmods+1)+'-'+str(j+1)+'\n Out 1'
                            btn['text'] = text
                            btn.bind("<Button-1>", self.ftest)
                            buttons[nmods].append(btn)
                        elif k == 2:
                            text = str(nmods+1)+'-'+str(j+1)+'\n Out 2'
                            btn['text'] = text
                            btn.bind("<Button-1>", self.ftest)
                            buttons[nmods].append(btn)
                nmods = nmods + 1

        #Place labels above the sum buttons in the grid.
        labels = []
        nlabels = 0
        for i in range(0,nslots):
            if i == 0 or i == 1 or i == 3 or i == 4 or i == 6 or i == 7 or i == 8 or i == 10 or i == 11:
                label = Label(RR2_lower_TDC_disc_frames[i], text='TDC Discriminator '+str(nlabels+1))
                labels.append(label)
                labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                    
        #Place the sum buttons in the grid.
        nmods = 0
        for i in range(0,nslots):
            if i == 0 or i == 1 or i == 3 or i == 4 or i == 6 or i == 7 or i == 8 or i == 10 or i == 11:
                for j in range(1,nrows+1):
                    for k in range(0,ncols):
                        buttons[nmods][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')
                nmods = nmods +1

    def RR3_splitter_connections(self, event):
        win = Tk()
        win.wm_title("RR3 Splitter Panel Connections")
        win.geometry("800x600")

        RR3_splitter_frames = []                      #List to hold the amplifier frames.
        buttons = []                                  #List to hold the amplifier buttons.
        nslots = 9                                    #Number of amplifiers in crate.
        nrows = 3                                     #Number of rows of channels in an amp.
        ncols = 16                                    #Number of columns of channels in an amp.
        nchannels = nrows * ncols                     #Number of channels in one amp.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the amps.

        canvas = Canvas(win, borderwidth=0)
        #Make a frame to go on the canvas that contains the other frames so the scroll bar works for all other frames.
        container_frame = Frame(canvas, relief=RAISED, borderwidth=1)
        container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        vsb = Scrollbar(win, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)

        vsb.pack(side="right", fill="y")
        canvas.pack(side="top", fill="both", expand=True)
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        canvas.create_window((0,0), window=container_frame, anchor="nw")
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Expands container frame but breaks scroll bar.
        container_frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure(canvas))

        #Create 9 different frames to hold the amp channel buttons.
        for i in range (0,nslots):
            frame = Frame(container_frame, relief=RAISED, borderwidth=2)
            frame.pack(side='top', fill=BOTH, expand=True)
            RR3_splitter_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            for j in range(0,ncols):
                RR3_splitter_frames[i].columnconfigure(j, pad=3, weight=1)

            for j in range(0,nrows+1):
                RR3_splitter_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the amp buttons to their functions.
        for i in range(0,nslots):
            buttons.append([])#Make this list 2D to hold the buttons of each of the 9 frames separately.
            for j in range(0,nrows):
                for k in range(0,ncols):

                    btn = Button(RR3_splitter_frames[i], width=3, height=2, font='Helvetica 8')

                    #Give names and function binds to the three different columns.
                    if j == 0:
                        text = 'SP'+str(18-i)+'-'+str(k+1)+'\n Out 1'
                        btn['text'] = text
                        #btn['font'] = 'Helvetica 20'
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)
                    elif j == 1:
                        text = 'SP'+str(18-i)+'-'+str(k+1)+'\n In'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)
                    elif j == 2:
                        text = 'SP'+str(18-i)+'-'+str(k+1)+'\n Out 2'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)

        labels = []
        #Place labels above the buttons in the grid.
        for i in range(0,nslots):
            label = Label(RR3_splitter_frames[i], text='Splitter Panel '+str(18-i))
            labels.append(label)
            labels[i].grid(row=0,columnspan=ncols)

        #Place the buttons in the grid.
        for i in range(0,nslots):
            for j in range(1,nrows+1):
                for k in range(0,ncols):
                    buttons[i][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')


        def onFrameConfigure(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            canvas.configure(scrollregion=canvas.bbox("all"))

    def RR1_splitter_connections(self, event):
        win = Tk()
        win.wm_title("RR1 Splitter Panel Connections")
        win.geometry("800x600")

        RR1_splitter_frames = []                      #List to hold the amplifier frames.
        buttons = []                                  #List to hold the amplifier buttons.
        nslots = 9                                    #Number of amplifiers in crate.
        nrows = 3                                     #Number of rows of channels in an amp.
        ncols = 16                                    #Number of columns of channels in an amp.
        nchannels = nrows * ncols                     #Number of channels in one amp.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the amps.

        canvas = Canvas(win, borderwidth=0)
        #Make a frame to go on the canvas that contains the other frames so the scroll bar works for all other frames.
        container_frame = Frame(canvas, relief=RAISED, borderwidth=1)
        container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        vsb = Scrollbar(win, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)

        vsb.pack(side="right", fill="y")
        canvas.pack(side="top", fill="both", expand=True)
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        canvas.create_window((0,0), window=container_frame, anchor="nw")
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Expands container frame but breaks scroll bar.
        container_frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure(canvas))

        #Create 9 different frames to hold the amp channel buttons.
        for i in range (0,nslots):
            frame = Frame(container_frame, relief=RAISED, borderwidth=2)
            frame.pack(side='top', fill=BOTH, expand=True)
            RR1_splitter_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            for j in range(0,ncols):
                RR1_splitter_frames[i].columnconfigure(j, pad=3, weight=1)

            for j in range(0,nrows+1):
                RR1_splitter_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the amp buttons to their functions.
        for i in range(0,nslots):
            buttons.append([])#Make this list 2D to hold the buttons of each of the 9 frames separately.
            for j in range(0,nrows):
                for k in range(0,ncols):

                    btn = Button(RR1_splitter_frames[i], width=3, height=2, font='Helvetica 8')

                    #Give names and function binds to the three different columns.
                    if j == 0:
                        text = 'SP'+str(nslots-i)+'-'+str(k+1)+'\n Out 1'
                        btn['text'] = text
                        #btn['font'] = 'Helvetica 20'
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)
                    elif j == 1:
                        text = 'SP'+str(nslots-i)+'-'+str(k+1)+'\n In'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)
                    elif j == 2:
                        text = 'SP'+str(nslots-i)+'-'+str(k+1)+'\n Out 2'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)

        labels = []
        #Place labels above the buttons in the grid.
        for i in range(0,nslots):
            label = Label(RR1_splitter_frames[i], text='Splitter Panel '+str(nslots-i))
            labels.append(label)
            labels[i].grid(row=0,columnspan=ncols)

        #Place the buttons in the grid.
        for i in range(0,nslots):
            for j in range(1,nrows+1):
                for k in range(0,ncols):
                    buttons[i][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')


        def onFrameConfigure(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            canvas.configure(scrollregion=canvas.bbox("all"))

    def RR2_adc_pp_connections(self, event):
        win = Tk()
        win.wm_title("RR2 fADC Patch Panel Connections")
        win.geometry("800x600")

        RR2_adc_pp_frames = []                        #List to hold the amplifier frames.
        buttons = []                                  #List to hold the amplifier buttons.
        nslots = 5                                    #Number of amplifiers in crate.
        nrows = 4                                     #Number of rows of channels in an amp.
        ncols = 16                                    #Number of columns of channels in an amp.
        nchannels = nrows * ncols                     #Number of channels in one amp.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the amps.

        canvas = Canvas(win, borderwidth=0)
        #Make a frame to go on the canvas that contains the other frames so the scroll bar works for all other frames.
        container_frame = Frame(canvas, relief=RAISED, borderwidth=1)
        container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        vsb = Scrollbar(win, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)

        vsb.pack(side="right", fill="y")
        canvas.pack(side="top", fill="both", expand=True)
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        canvas.create_window((0,0), window=container_frame, anchor="nw")
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Expands container frame but breaks scroll bar.
        container_frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure2(canvas))

        #Create 9 different frames to hold the amp channel buttons.
        for i in range (0,nslots):
            frame = Frame(container_frame, relief=RAISED, borderwidth=2)
            frame.pack(side='top', fill=BOTH, expand=True)
            RR2_adc_pp_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            for j in range(0,ncols):
                RR2_adc_pp_frames[i].columnconfigure(j, pad=3, weight=1)

            for j in range(0,nrows+1):
                RR2_adc_pp_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the amp buttons to their functions.
        for i in range(0,nslots):
            buttons.append([])#Make this list 2D to hold the buttons of each of the 9 frames separately.
            for j in range(0,nrows):
                for k in range(0,ncols):

                    btn = Button(RR2_adc_pp_frames[i], width=3, height=2, font='Helvetica 8')

                    #Give names and function binds to the three different columns.
                    if j == 0:
                        text = str(nslots-i)+'A-'+str(k+1)
                        btn['text'] = text
                        #btn['font'] = 'Helvetica 20'
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)
                    elif j == 1:
                        text = str(nslots-i)+'B-'+str(k+1)
                        btn['text'] = text
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)
                    elif j == 2:
                        text = str(nslots-i)+'C-'+str(k+1)
                        btn['text'] = text
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)
                    elif j == 3:
                        text = str(nslots-i)+'D-'+str(k+1)
                        btn['text'] = text
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)

        labels = []
        #Place labels above the buttons in the grid.
        for i in range(0,nslots):
            label = Label(RR2_adc_pp_frames[i], text='fADC Patch Panel '+str(nslots-i))
            labels.append(label)
            labels[i].grid(row=0,columnspan=ncols)

        #Place the buttons in the grid.
        for i in range(0,nslots):
            for j in range(1,nrows+1):
                for k in range(0,ncols):
                    buttons[i][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')


        def onFrameConfigure2(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            canvas.configure(scrollregion=canvas.bbox("all"))

    def RR2_lower_tdc_pp_connections(self, event):
        win = Tk()
        win.wm_title("RR2 Lower TDC Patch Panel Connections")
        win.geometry("800x600")

        RR2_lower_tdc_pp_frames = []                  #List to hold the amplifier frames.
        buttons = []                                  #List to hold the amplifier buttons.
        nslots = 3                                    #Number of amplifiers in crate.
        nrows = 4                                     #Number of rows of channels in an amp.
        ncols = 16                                    #Number of columns of channels in an amp.
        nchannels = nrows * ncols                     #Number of channels in one amp.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the amps.

        canvas = Canvas(win, borderwidth=0)
        #Make a frame to go on the canvas that contains the other frames so the scroll bar works for all other frames.
        container_frame = Frame(canvas, relief=RAISED, borderwidth=1)
        container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        vsb = Scrollbar(win, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)

        vsb.pack(side="right", fill="y")
        canvas.pack(side="top", fill="both", expand=True)
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        canvas.create_window((0,0), window=container_frame, anchor="nw")
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Expands container frame but breaks scroll bar.
        container_frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure2(canvas))

        #Create 9 different frames to hold the amp channel buttons.
        for i in range (0,nslots):
            frame = Frame(container_frame, relief=RAISED, borderwidth=2)
            frame.pack(side='top', fill=BOTH, expand=True)
            RR2_lower_tdc_pp_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            for j in range(0,ncols):
                RR2_lower_tdc_pp_frames[i].columnconfigure(j, pad=3, weight=1)

            for j in range(0,nrows+1):
                RR2_lower_tdc_pp_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the amp buttons to their functions.
        for i in range(0,nslots):
            buttons.append([])#Make this list 2D to hold the buttons of each of the 9 frames separately.
            for j in range(0,nrows):
                for k in range(0,ncols):

                    btn = Button(RR2_lower_tdc_pp_frames[i], width=3, height=2, font='Helvetica 8')

                    #Give names and function binds to the three different columns.
                    if j == 0:
                        text = str(nslots-i)+'A-'+str(k+1)
                        btn['text'] = text
                        #btn['font'] = 'Helvetica 20'
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)
                    elif j == 1:
                        text = str(nslots-i)+'B-'+str(k+1)
                        btn['text'] = text
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)
                    elif j == 2:
                        text = str(nslots-i)+'C-'+str(k+1)
                        btn['text'] = text
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)
                    elif j == 3:
                        text = str(nslots-i)+'D-'+str(k+1)
                        btn['text'] = text
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)

        labels = []
        #Place labels above the buttons in the grid.
        for i in range(0,nslots):
            label = Label(RR2_lower_tdc_pp_frames[i], text='fADC Patch Panel '+str(nslots-i))
            labels.append(label)
            labels[i].grid(row=0,columnspan=ncols)

        #Place the buttons in the grid.
        for i in range(0,nslots):
            for j in range(1,nrows+1):
                for k in range(0,ncols):
                    buttons[i][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')


        def onFrameConfigure2(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            canvas.configure(scrollregion=canvas.bbox("all"))

    def RR2_upper_tdc_pp_connections(self, event):
        win = Tk()
        win.wm_title("RR2 Upper TDC Patch Panel Connections")
        win.geometry("800x600")

        RR2_upper_tdc_pp_frames = []                  #List to hold the amplifier frames.
        buttons = []                                  #List to hold the amplifier buttons.
        nslots = 2                                    #Number of amplifiers in crate.
        nrows = 4                                     #Number of rows of channels in an amp.
        ncols = 16                                    #Number of columns of channels in an amp.
        nchannels = nrows * ncols                     #Number of channels in one amp.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the amps.

        canvas = Canvas(win, borderwidth=0)
        #Make a frame to go on the canvas that contains the other frames so the scroll bar works for all other frames.
        container_frame = Frame(canvas, relief=RAISED, borderwidth=1)
        container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        vsb = Scrollbar(win, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)

        vsb.pack(side="right", fill="y")
        canvas.pack(side="top", fill="both", expand=True)
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        canvas.create_window((0,0), window=container_frame, anchor="nw")
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Expands container frame but breaks scroll bar.
        container_frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure2(canvas))

        #Create 9 different frames to hold the amp channel buttons.
        for i in range (0,nslots):
            frame = Frame(container_frame, relief=RAISED, borderwidth=2)
            frame.pack(side='top', fill=BOTH, expand=True)
            RR2_upper_tdc_pp_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            for j in range(0,ncols):
                RR2_upper_tdc_pp_frames[i].columnconfigure(j, pad=3, weight=1)

            for j in range(0,nrows+1):
                RR2_upper_tdc_pp_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the amp buttons to their functions.
        for i in range(0,nslots):
            buttons.append([])#Make this list 2D to hold the buttons of each of the 9 frames separately.
            for j in range(0,nrows):
                for k in range(0,ncols):

                    btn = Button(RR2_upper_tdc_pp_frames[i], width=3, height=2, font='Helvetica 8')

                    #Give names and function binds to the three different columns.
                    if j == 0:
                        text = str(5-i)+'A-'+str(k+1)
                        btn['text'] = text
                        #btn['font'] = 'Helvetica 20'
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)
                    elif j == 1:
                        text = str(5-i)+'B-'+str(k+1)
                        btn['text'] = text
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)
                    elif j == 2:
                        text = str(5-i)+'C-'+str(k+1)
                        btn['text'] = text
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)
                    elif j == 3:
                        text = str(5-i)+'D-'+str(k+1)
                        btn['text'] = text
                        btn.bind("<Button-1>", self.ftest)
                        buttons[i].append(btn)

        labels = []
        #Place labels above the buttons in the grid.
        for i in range(0,nslots):
            label = Label(RR2_upper_tdc_pp_frames[i], text='fADC Patch Panel '+str(5-i))
            labels.append(label)
            labels[i].grid(row=0,columnspan=ncols)

        #Place the buttons in the grid.
        for i in range(0,nslots):
            for j in range(1,nrows+1):
                for k in range(0,ncols):
                    buttons[i][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')


        def onFrameConfigure2(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            canvas.configure(scrollregion=canvas.bbox("all"))

root = Tk()
root.geometry("1200x800")
my_gui = MyFirstGUI(root)
root.mainloop()
