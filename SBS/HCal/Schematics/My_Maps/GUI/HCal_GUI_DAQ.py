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

class MyFirstGUI:
    def __init__(self, primary):
        self.primary = primary
        primary.title("Front-End Electronics")

        self.style = Style()
        self.style.theme_use("default")

        #Create three frames to hold the contents of each DAQ rack.
        RR4_frame = Frame(primary, relief=RAISED, borderwidth=2)
        RR4_frame.pack(side='left', fill=BOTH, expand=True)

        RR5_frame = Frame(primary, relief=RAISED, borderwidth=2)
        RR5_frame.pack(side='left', fill=BOTH, expand=True)

        #Create buttons and labels for electronics in RR4.
        RR4_label = Label(RR4_frame, text='RR4', font='Helvetica 16 bold')
        RR4_label.pack(side="top")

        RR4_tdc_pp_btn = Button(RR4_frame, text='TDC Patch Panels', width=10, height=6, font='Helvetica 14 bold')
        RR4_tdc_pp_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR4_tdc_pp_btn.bind("<Button-1>", self.RR4_tdc_pp_connections)

        RR4_tdc_disc_btn = Button(RR4_frame, text='TDC Discriminators', width=10, height=3, font='Helvetica 14 bold')
        RR4_tdc_disc_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR4_tdc_disc_btn.bind("<Button-1>", self.RR4_tdc_disc_connections)

        RR4_adc_pp_btn = Button(RR4_frame, text='ADC Patch Panels', width=10, height=6, font='Helvetica 14 bold')
        RR4_adc_pp_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR4_adc_pp_btn.bind("<Button-1>", self.RR4_adc_pp_connections)

        #Create buttons and labels for electronics in RR5.
        RR5_label = Label(RR5_frame, text='RR5', font='Helvetica 16 bold')
        RR5_label.pack(side="top")

        RR5_nim_btn = Button(RR5_frame, text='NIM Electronics', width=10, height=3, font='Helvetica 14 bold')
        RR5_nim_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        #RR5_nim_btn.bind("<Button-1>", self.RR5_tdc_pp_connections)

        RR5_vxs2_btn = Button(RR5_frame, text='VXS 2 (intelsbshcal2)', width=10, height=6, font='Helvetica 14 bold')
        RR5_vxs2_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR5_vxs2_btn.bind("<Button-1>", self.RR5_vxs2_connections)

        RR5_vxs1_btn = Button(RR5_frame, text='VXS 1 (intelsbshcal1)', width=10, height=6, font='Helvetica 14 bold')
        RR5_vxs1_btn.pack(side="top", padx=5, pady=5, expand=True, fill=BOTH)
        RR5_vxs1_btn.bind("<Button-1>", self.RR5_vxs1_connections)

    def ftest(self, event):
        print('Test successful!', event.widget.cget('text'))

    def RR4_tdc_pp_connections(self, event):
        win = Tk()
        win.wm_title("RR2 Upper TDC Patch Panel Connections")
        win.geometry("800x600")

        RR4_tdc_pp_frames = []                  #List to hold the amplifier frames.
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
            RR4_tdc_pp_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            for j in range(0,ncols):
                RR4_tdc_pp_frames[i].columnconfigure(j, pad=3, weight=1)

            for j in range(0,nrows+1):
                RR4_tdc_pp_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the amp buttons to their functions.
        for i in range(0,nslots):
            buttons.append([])#Make this list 2D to hold the buttons of each of the 9 frames separately.
            for j in range(0,nrows):
                for k in range(0,ncols):

                    btn = Button(RR4_tdc_pp_frames[i], width=3, height=2, font='Helvetica 8')

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
            label = Label(RR4_tdc_pp_frames[i], text='TDC Patch Panel '+str(5-i))
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

    def RR4_tdc_disc_connections(self, event):
        win = Tk()
        win.wm_title("RR2 Upper TDC Discriminator Connections")
        win.geometry("1400x525")
        RR4_tdc_disc_frames = []                #List to hold the frames.
        buttons = []                                  #List to hold the  buttons.
        nslots = 25                                   #Number of slots in crate.
        ndisc = 18                                     #Number of summing modules.
        nrows = 8                                    #Number of rows of channels in a module.
        ncols = 2                                     #Number of columns of channels in a module.
        nchannels = nrows * ncols                     #Number of channels in one module.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the modules.

        canvas = Canvas(win, borderwidth=0)
        #Make a frame to go on the canvas that contains the other frames so the scroll bar works for all other frames.
        container_frame = Frame(canvas, relief=RAISED, borderwidth=1)
        container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        vsb = Scrollbar(win, orient="horizontal", command=canvas.xview)
        canvas.configure(yscrollcommand=vsb.set)

        vsb.pack(side="bottom", fill="x")
        canvas.pack(side="top", fill="both", expand=True)
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        canvas.create_window((0,0), window=container_frame, anchor="nw")
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Expands container frame but breaks scroll bar.
        container_frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure2(canvas))

        #Create different frames to hold the channel buttons.
        for i in range (0,nslots):
            frame = Frame(container_frame, relief=RAISED, borderwidth=2)
            frame.pack(side='left', fill=BOTH, expand=True)
            RR4_tdc_disc_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            if i == 1 or i == 2 or i == 3 or i == 5 or i == 6 or i == 7 or i == 9 or i == 10 or i == 11 or i == 13 or i == 14 or i == 15 or i == 17 or i == 18 or i == 19 or i == 21 or i == 22 or i == 23 :
                for j in range(0,ncols):
                    RR4_tdc_disc_frames[i].columnconfigure(j, pad=3, weight=1)

                for j in range(0,nrows+3):
                    RR4_tdc_disc_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the sum buttons to their functions.
        nmods = 0
        for i in range(0,nslots):
            if i == 1 or i == 2 or i == 3 or i == 5 or i == 6 or i == 7 or i == 9 or i == 10 or i == 11 or i == 13 or i == 14 or i == 15 or i == 17 or i == 18 or i == 19 or i == 21 or i == 22 or i == 23 :
                buttons.append([])#Make this list 2D to hold the buttons of each of the frames separately.
                for j in range(0,nrows):
                    for k in range(0,ncols):
                        btn = Button(RR4_tdc_disc_frames[i], width=1, height=1, font='Helvetica 8')

                        #Give names and function binds to the three different columns.
                        if k == 0:
                            text = str(18-nmods)+'-'+str(j*2+1)+'\n In'
                            btn['text'] = text
                            btn.bind("<Button-1>", self.ftest)
                            buttons[nmods].append(btn)
                        elif k == 1:
                            text = str(18-nmods)+'-'+str(j*2+2)+'\n In'
                            btn['text'] = text
                            btn.bind("<Button-1>", self.ftest)
                            buttons[nmods].append(btn)
                #Make a button for each discriminator output.
                btn1 = Button(RR4_tdc_disc_frames[i], width=1, height=8, font='Helvetica 8')
                btn1['text'] = 'Ribbon \n Out 1'
                btn1.bind("<Button-1>", self.ftest)
                buttons[nmods].append(btn1)
                btn2 = Button(RR4_tdc_disc_frames[i], width=1, height=8, font='Helvetica 8')
                btn2['text'] = 'Ribbon \n Out 2'
                btn2.bind("<Button-1>", self.ftest)
                buttons[nmods].append(btn2)

                nmods = nmods + 1

        #Place labels above the sum buttons in the grid.
        labels = []
        nlabels = 0
        nmod = 0
        for i in range(0,nslots):
            if i == 1 or i == 2 or i == 3 or i == 5 or i == 6 or i == 7 or i == 9 or i == 10 or i == 11 or i == 13 or i == 14 or i == 15 or i == 17 or i == 18 or i == 19 or i == 21 or i == 22 or i == 23 :
                label = Label(RR4_tdc_disc_frames[i], text='TDC Discriminator '+str(ndisc-nmod))
                labels.append(label)
                labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                nmod = nmod + 1
            else:
                label = Label(RR4_tdc_disc_frames[i], text='Empty')
                labels.append(label)
                labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                    
        #Place the buttons in the grid.
        nmods = 0
        for i in range(0,nslots):
            if i == 1 or i == 2 or i == 3 or i == 5 or i == 6 or i == 7 or i == 9 or i == 10 or i == 11 or i == 13 or i == 14 or i == 15 or i == 17 or i == 18 or i == 19 or i == 21 or i == 22 or i == 23 :
                for j in range(1,5):
                    for k in range(0,ncols):
                        buttons[nmods][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')

                for j in range(7,11):
                    for k in range(0,ncols):
                        buttons[nmods][(j-1-2)*ncols+k].grid(row=j, column=k, sticky='NSEW')
                        
                #Add the two output buttons.
                buttons[nmods][16].grid(row=5, column=0, rowspan=1, sticky='NSEW')
                buttons[nmods][17].grid(row=6, column=0, rowspan=1, sticky='NSEW')
                nmods = nmods +1


        def onFrameConfigure2(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            canvas.configure(scrollregion=canvas.bbox("all"))

    def RR4_adc_pp_connections(self, event):
        win = Tk()
        win.wm_title("RR2 Upper TDC Patch Panel Connections")
        win.geometry("800x600")

        RR4_adc_pp_frames = []                  #List to hold the amplifier frames.
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
            RR4_adc_pp_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            for j in range(0,ncols):
                RR4_adc_pp_frames[i].columnconfigure(j, pad=3, weight=1)

            for j in range(0,nrows+1):
                RR4_adc_pp_frames[i].rowconfigure(j, pad=3, weight=1)

        #Create and bind the amp buttons to their functions.
        for i in range(0,nslots):
            buttons.append([])#Make this list 2D to hold the buttons of each of the 9 frames separately.
            for j in range(0,nrows):
                for k in range(0,ncols):

                    btn = Button(RR4_adc_pp_frames[i], width=3, height=2, font='Helvetica 8')

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
            label = Label(RR4_adc_pp_frames[i], text='fADC Patch Panel '+str(5-i))
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

    def RR5_vxs2_connections(self, event):
        win = Tk()
        win.wm_title("RR5 Upper VXS Crate Connections")
        win.geometry("1400x705")
        RR5_vxs2_frames = []                #List to hold the frames.
        buttons = []                                  #List to hold the  buttons.
        ti_buttons = []
        f1_buttons = []
        nslots = 21                                   #Number of slots in crate.
        ndisc = 2                                     #Number of summing modules.
        nrows = 16                                    #Number of rows of channels in a module.
        ncols = 1                                     #Number of columns of channels in a module.
        ti_rows = 18
        ti_cols = 1
        nf1 = 5
        f1_rows = 3
        f1_cols = 2
        nchannels = nrows * ncols                     #Number of channels in one module.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the modules.

        canvas = Canvas(win, borderwidth=0)
        #Make a frame to go on the canvas that contains the other frames so the scroll bar works for all other frames.
        container_frame = Frame(canvas, relief=RAISED, borderwidth=1)
        container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        vsb = Scrollbar(win, orient="horizontal", command=canvas.xview)
        canvas.configure(yscrollcommand=vsb.set)

        vsb.pack(side="bottom", fill="x")
        canvas.pack(side="top", fill="both", expand=True)
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        canvas.create_window((0,0), window=container_frame, anchor="nw")
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Expands container frame but breaks scroll bar.
        container_frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure2(canvas))

        #Create different frames to hold the channel buttons.
        for i in range (0,nslots):
            frame = Frame(container_frame, relief=RAISED, borderwidth=2)
            frame.pack(side='left', fill=BOTH, expand=True)
            RR5_vxs2_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            if i == 18 or i == 19:
                for j in range(0,ncols):
                    RR5_vxs2_frames[i].columnconfigure(j, pad=3, weight=1)

                for j in range(0,nrows+1):
                    RR5_vxs2_frames[i].rowconfigure(j, pad=3, weight=1)

        #Define standalone module rows and cols.
        RR5_vxs2_frames[0].columnconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[0].rowconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[0].rowconfigure(1, pad=3, weight=1)

        RR5_vxs2_frames[1].columnconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[1].rowconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[1].rowconfigure(1, pad=3, weight=1)

        RR5_vxs2_frames[10].columnconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[10].rowconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[10].rowconfigure(1, pad=3, weight=1)

        RR5_vxs2_frames[11].columnconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[11].rowconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[11].rowconfigure(1, pad=3, weight=1)

        RR5_vxs2_frames[14].columnconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[14].rowconfigure(0, pad=3, weight=1)
        RR5_vxs2_frames[14].rowconfigure(1, pad=3, weight=1)

        for i in range(0,nslots):
            if i == 2 or i == 8  or i == 9  or i == 12  or i == 13  or i == 15 or i == 16 or i == 17:
                RR5_vxs2_frames[i].columnconfigure(0, pad=3, weight=1)
                RR5_vxs2_frames[i].rowconfigure(0, pad=3, weight=1)
                RR5_vxs2_frames[i].rowconfigure(1, pad=3, weight=1)

        #Define TI button grid.
        for i in range(0,ti_rows+5):
            RR5_vxs2_frames[20].rowconfigure(i, pad=3, weight=1)

        for i in range(0,ti_cols):
            RR5_vxs2_frames[20].columnconfigure(i, pad=3, weight=1)

        #Define F1 button grid.
        for i in range(0,nslots):
            if i == 3 or i == 4 or i == 5 or i == 6 or i == 7:
                for j in range(0,f1_rows+1):
                    RR5_vxs2_frames[i].rowconfigure(j, pad=3, weight=1)

                for j in range(0,f1_cols):
                    RR5_vxs2_frames[i].columnconfigure(j, pad=3, weight=1)

        #Create and bind the fADC buttons to their functions.
        nmods = 0
        for i in range(0,nslots):
            if  i == 18 or i == 19:
                buttons.append([])#Make this list 2D to hold the buttons of each of the frames separately.
                for j in range(0,nrows):
                    for k in range(0,ncols):
                        btn = Button(RR5_vxs2_frames[i], width=1, height=1, font='Helvetica 8')

                        #Give names and function binds to the three different columns.
                        text = str(nmods+1)+'-'+str(15-j)+'\n In'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.ftest)
                        buttons[nmods].append(btn)
                nmods = nmods + 1

        #Create and bind the F1 buttons to their functions.
        nmods = 0
        for i in range(0,nslots):
            if i == 3 or i == 4 or i == 5 or i == 6 or i == 7:
                f1_buttons.append([])#Make this list 2D to hold the buttons of each of the frames separately.
                for j in range(0,f1_rows):
                    for k in range(0,f1_cols):
                        f1_btn = Button(RR5_vxs2_frames[i], width=1, height=10, font='Helvetica 8')

                        #Give names and function binds to the three different columns.
                        if j == 0 and k == 0:
                            text = 'Chs In\n'+str(17)+'-'+str(32)
                        if j == 0 and k == 1:
                            text = 'Chs In\n'+str(1)+'-'+str(16)
                        if j == 1 and k == 0:
                            text = 'F1 SD\n Cable'
                            f1_btn['height'] = 1
                        if j == 1 and k == 1:
                            text = ''
                            f1_btn['height'] = 1
                            f1_btn['state'] = 'disabled'
                        if j == 2 and k == 0:
                            text = 'Chs In\n'+str(49)+'-'+str(64)
                        if j == 2 and k == 1:
                            text = 'Chs In\n'+str(33)+'-'+str(48)
                        f1_btn['text'] = text
                        f1_btn.bind("<Button-1>", self.ftest)
                        f1_buttons[nmods].append(f1_btn)
                nmods = nmods + 1

        #Create a standalone buttons.
        roc_btn = Button(RR5_vxs2_frames[0], text='ROC 17', width=4, height=49, font='Helvetica 8')
        roc_btn.grid(row=1,columnspan=ncols)
        roc_btn.bind("<Button-1>", self.ftest)

        scaler_btn = Button(RR5_vxs2_frames[1], text='F1TDC SD', width=4, height=49, font='Helvetica 8')
        scaler_btn.grid(row=1,columnspan=ncols)
        scaler_btn.bind("<Button-1>", self.ftest)

        vtp_btn = Button(RR5_vxs2_frames[10], text='VTP', width=4, height=49, font='Helvetica 8')
        vtp_btn.grid(row=1,columnspan=ncols)
        vtp_btn.bind("<Button-1>", self.ftest)

        fadc_sd_btn = Button(RR5_vxs2_frames[11], text='fADC SD', width=4, height=49, font='Helvetica 8')
        fadc_sd_btn.grid(row=1,columnspan=ncols)
        fadc_sd_btn.bind("<Button-1>", self.ftest)

        dvcs_pulser_btn = Button(RR5_vxs2_frames[14], text='DVCS\n Pulser', width=4, height=49, font='Helvetica 8')
        dvcs_pulser_btn.grid(row=1,columnspan=ncols)
        dvcs_pulser_btn.bind("<Button-1>", self.ftest)

        #Create buttons for empty slots.
        empty_btns = []
        for i in range(0,nslots):
            if i == 2 or i == 8  or i == 9  or i == 12  or i == 13  or i == 15 or i == 16 or i == 17:
                empty_btn = Button(RR5_vxs2_frames[i], text='Empty', width=4, height=49, font='Helvetica 8')
                empty_btn.grid(row=1,columnspan=ncols)
                empty_btn.bind("<Button-1>", self.ftest)
                empty_btn['state'] = 'disabled'
                empty_btns.append(empty_btn)

        #Make TI buttons.
        ti_btn_names = ['Optical\n Fiber','','','','','TS#6\n In','TS#5\n In','TS#4\n In','CLK\n In','TS#3\n In','TS#2\n In','TS#1\n In','TRG\n In','INH\n In','OT#2\n Out','OT#1\n Out','O#4\n Out','O#3\n Out','O#2\n Out','O#1\n Out','TRG\n Out','Busy\n Out']

        for i in range(0,ti_rows+4):
            for j in range(0,ti_cols):
                ti_btn = Button(RR5_vxs2_frames[20], width=1, height=1, font='Helvetica 8')
                ti_btn['text'] = ti_btn_names[i]
                ti_btn.bind("<Button-1>", self.ftest)

                if i == 1 or i == 2 or i == 3 or i == 4:
                    ti_btn['state'] = 'disabled'

                ti_buttons.append(ti_btn)

        #Place labels above the fADC buttons in the grid.
        labels = []
        nlabels = 0
        nmod = 0
        for i in range(0,nslots):
            if i == 18 or i == 19:
                label = Label(RR5_vxs2_frames[i], text='fADC250 '+str(i-1))
                labels.append(label)
                labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                nmod = nmod + 1

        #Place labels above the F1 buttons in the grid.
        f1_labels = []
        nlabels = 0
        nmod = 0
        for i in range(0,nslots):
           if i == 3 or i == 4 or i == 5 or i == 6 or i == 7:
                f1_label = Label(RR5_vxs2_frames[i], text='F1 TDC '+str(i-2))
                f1_labels.append(f1_label)
                f1_labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                nmod = nmod + 1
              
        #Place labels above the empty buttons in the grid.
        empty_labels = []
        nlabels = 0
        nmod = 0
        for i in range(0,nslots):
            if i == 2 or i == 8  or i == 9  or i == 12  or i == 13  or i == 15 or i == 16 or i == 17:
                empty_label = Label(RR5_vxs2_frames[i], text='Empty')
                empty_labels.append(empty_label)
                empty_labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                nmod = nmod + 1
                

        #Create standalone labels.
        roc_label = Label(RR5_vxs2_frames[0], text='ROC17')
        roc_label.grid(row=0,columnspan=ncols)

        scaler_label = Label(RR5_vxs2_frames[1], text='F1TDC SD')
        scaler_label.grid(row=0,columnspan=ncols)

        vtp_label = Label(RR5_vxs2_frames[10], text='VTP')
        vtp_label.grid(row=0,columnspan=ncols)

        fadc_sd_label = Label(RR5_vxs2_frames[11], text='fADC SD')
        fadc_sd_label.grid(row=0,columnspan=ncols)

        ti_label = Label(RR5_vxs2_frames[20], text='TI')
        ti_label.grid(row=0,columnspan=ncols)

        dvcs_pulser_label = Label(RR5_vxs2_frames[14], text='DVCS Pulser')
        dvcs_pulser_label.grid(row=0,columnspan=ncols)

        #Place the fADC buttons in the grid.
        nmods = 0
        for i in range(0,nslots):
            if i == 18 or i == 19:
                for j in range(1,nrows+1):
                    for k in range(0,ncols):
                        buttons[nmods][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')
                nmods = nmods +1

        #Place the TI buttons in the grid.
        for i in range(1,ti_rows+5):
            for j in range(0,ti_cols):
                ti_buttons[i-1].grid(row=i, column=j, sticky='NSEW')

        #Place the f1 buttons in the grid.
        nmods = 0
        for i in range(0,nslots):
           if i == 3 or i == 4 or i == 5 or i == 6 or i == 7:
                for j in range(1,f1_rows+1):
                    for k in range(0,f1_cols):
                        f1_buttons[nmods][(j-1)*f1_cols+k].grid(row=j, column=k, sticky='NSEW')
                nmods = nmods +1


        def onFrameConfigure2(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            canvas.configure(scrollregion=canvas.bbox("all"))

    def RR5_vxs1_connections(self, event):
        win = Tk()
        win.wm_title("RR5 Lower VXS Crate Connections")
        win.geometry("1400x705")
        RR5_vxs1_frames = []                #List to hold the frames.
        buttons = []                                  #List to hold the  buttons.
        ti_buttons = []
        nslots = 21                                   #Number of slots in crate.
        ndisc = 16                                     #Number of summing modules.
        nrows = 16                                    #Number of rows of channels in a module.
        ncols = 1                                     #Number of columns of channels in a module.
        ti_rows = 18
        ti_cols = 1
        nchannels = nrows * ncols                     #Number of channels in one module.
        nall_channels = nslots * nrows * ncols        #Number of channels in all the modules.

        canvas = Canvas(win, borderwidth=0)
        #Make a frame to go on the canvas that contains the other frames so the scroll bar works for all other frames.
        container_frame = Frame(canvas, relief=RAISED, borderwidth=1)
        container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        vsb = Scrollbar(win, orient="horizontal", command=canvas.xview)
        canvas.configure(yscrollcommand=vsb.set)

        vsb.pack(side="bottom", fill="x")
        canvas.pack(side="top", fill="both", expand=True)
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Cannot for life of me get this to expand to entire width of popup. Currently just changing geometry to make it fit.
        canvas.create_window((0,0), window=container_frame, anchor="nw")
        #container_frame.pack(side='top', fill=BOTH, expand=True)#Expands container frame but breaks scroll bar.
        container_frame.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure2(canvas))

        #Create different frames to hold the channel buttons.
        for i in range (0,nslots):
            frame = Frame(container_frame, relief=RAISED, borderwidth=2)
            frame.pack(side='left', fill=BOTH, expand=True)
            RR5_vxs1_frames.append(frame)

        #Define the rows and cols for the button grid.
        for i in range(0,nslots):
            if i == 2 or i == 3 or i == 4 or i == 5 or i == 6 or i == 7 or i == 8 or i == 9 or i == 12 or i == 13 or i == 14 or i == 15 or i == 16 or i == 17 or i == 18 or i == 19:
                for j in range(0,ncols):
                    RR5_vxs1_frames[i].columnconfigure(j, pad=3, weight=1)

                for j in range(0,nrows+1):
                    RR5_vxs1_frames[i].rowconfigure(j, pad=3, weight=1)
        #Define standalone module rows and cols.
        RR5_vxs1_frames[0].columnconfigure(0, pad=3, weight=1)
        RR5_vxs1_frames[0].rowconfigure(0, pad=3, weight=1)
        RR5_vxs1_frames[0].rowconfigure(1, pad=3, weight=1)

        RR5_vxs1_frames[1].columnconfigure(0, pad=3, weight=1)
        RR5_vxs1_frames[1].rowconfigure(0, pad=3, weight=1)
        RR5_vxs1_frames[1].rowconfigure(1, pad=3, weight=1)

        RR5_vxs1_frames[10].columnconfigure(0, pad=3, weight=1)
        RR5_vxs1_frames[10].rowconfigure(0, pad=3, weight=1)
        RR5_vxs1_frames[10].rowconfigure(1, pad=3, weight=1)

        RR5_vxs1_frames[11].columnconfigure(0, pad=3, weight=1)
        RR5_vxs1_frames[11].rowconfigure(0, pad=3, weight=1)
        RR5_vxs1_frames[11].rowconfigure(1, pad=3, weight=1)

        #Define TI button grid.
        for i in range(0,ti_rows+5):
            RR5_vxs1_frames[20].rowconfigure(i, pad=3, weight=1)

        for i in range(0,ti_cols):
            RR5_vxs1_frames[20].columnconfigure(i, pad=3, weight=1)

        #Create and bind the sum buttons to their functions.
        nmods = 0
        for i in range(0,nslots):
            if i == 2 or i == 3 or i == 4 or i == 5 or i == 6 or i == 7 or i == 8 or i == 9 or i == 12 or i == 13 or i == 14 or i == 15 or i == 16 or i == 17 or i == 18 or i == 19:
                buttons.append([])#Make this list 2D to hold the buttons of each of the frames separately.
                for j in range(0,nrows):
                    for k in range(0,ncols):
                        btn = Button(RR5_vxs1_frames[i], width=1, height=1, font='Helvetica 8')

                        #Give names and function binds to the three different columns.
                        text = str(nmods+1)+'-'+str(15-j)+'\n In'
                        btn['text'] = text
                        btn.bind("<Button-1>", self.ftest)
                        buttons[nmods].append(btn)
                nmods = nmods + 1
        #Create a standalone buttons.
        roc_btn = Button(RR5_vxs1_frames[0], text='ROC 16', width=4, height=49, font='Helvetica 8')
        roc_btn.grid(row=1,columnspan=ncols)
        roc_btn.bind("<Button-1>", self.ftest)

        scaler_btn = Button(RR5_vxs1_frames[1], text='Scaler', width=4, height=49, font='Helvetica 8')
        scaler_btn.grid(row=1,columnspan=ncols)
        scaler_btn.bind("<Button-1>", self.ftest)

        vtp_btn = Button(RR5_vxs1_frames[10], text='VTP', width=4, height=49, font='Helvetica 8')
        vtp_btn.grid(row=1,columnspan=ncols)
        vtp_btn.bind("<Button-1>", self.ftest)

        fadc_sd_btn = Button(RR5_vxs1_frames[11], text='fADC SD', width=4, height=49, font='Helvetica 8')
        fadc_sd_btn.grid(row=1,columnspan=ncols)
        fadc_sd_btn.bind("<Button-1>", self.ftest)

        #Make TI buttons.
        ti_btn_names = ['Optical\n Fiber','','','','','TS#6\n In','TS#5\n In','TS#4\n In','CLK\n In','TS#3\n In','TS#2\n In','TS#1\n In','TRG\n In','INH\n In','OT#2\n Out','OT#1\n Out','O#4\n Out','O#3\n Out','O#2\n Out','O#1\n Out','TRG\n Out','Busy\n Out']

        for i in range(0,ti_rows+4):
            for j in range(0,ti_cols):
                ti_btn = Button(RR5_vxs1_frames[20], width=1, height=1, font='Helvetica 8')
                ti_btn['text'] = ti_btn_names[i]
                ti_btn.bind("<Button-1>", self.ftest)

                if i == 1 or i == 2 or i == 3 or i == 4:
                    ti_btn['state'] = 'disabled'

                ti_buttons.append(ti_btn)

        #Place labels above the sum buttons in the grid.
        labels = []
        nlabels = 0
        nmod = 0
        for i in range(0,nslots):
            if i == 2 or i == 3 or i == 4 or i == 5 or i == 6 or i == 7 or i == 8 or i == 9 or i == 12 or i == 13 or i == 14 or i == 15 or i == 16 or i == 17 or i == 18 or i == 19:
                label = Label(RR5_vxs1_frames[i], text='fADC250 '+str(nmod+1))
                labels.append(label)
                labels[nlabels].grid(row=0,columnspan=ncols)
                nlabels = nlabels +1
                nmod = nmod + 1
                 
        #Create standalone labels.
        roc_label = Label(RR5_vxs1_frames[0], text='ROC16')
        roc_label.grid(row=0,columnspan=ncols)

        scaler_label = Label(RR5_vxs1_frames[1], text='Scaler')
        scaler_label.grid(row=0,columnspan=ncols)

        vtp_label = Label(RR5_vxs1_frames[10], text='VTP')
        vtp_label.grid(row=0,columnspan=ncols)

        fadc_sd_label = Label(RR5_vxs1_frames[11], text='fADC SD')
        fadc_sd_label.grid(row=0,columnspan=ncols)

        ti_label = Label(RR5_vxs1_frames[20], text='TI')
        ti_label.grid(row=0,columnspan=ncols)

        #Place the fADC buttons in the grid.
        nmods = 0
        for i in range(0,nslots):
            if i == 2 or i == 3 or i == 4 or i == 5 or i == 6 or i == 7 or i == 8 or i == 9 or i == 12 or i == 13 or i == 14 or i == 15 or i == 16 or i == 17 or i == 18 or i == 19:
                for j in range(1,nrows+1):
                    for k in range(0,ncols):
                        buttons[nmods][(j-1)*ncols+k].grid(row=j, column=k, sticky='NSEW')
                nmods = nmods +1

        #Place the TI buttons in the grid.
        for i in range(1,ti_rows+5):
            for j in range(0,ti_cols):
                ti_buttons[i-1].grid(row=i, column=j, sticky='NSEW')

        def onFrameConfigure2(canvas):
            '''Reset the scroll region to encompass the inner frame'''
            canvas.configure(scrollregion=canvas.bbox("all"))

root = Tk()
root.geometry("1200x800")
my_gui = MyFirstGUI(root)
root.mainloop()