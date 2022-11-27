from tkinter import Tk, Label, Button, StringVar, W, E, Frame, RIGHT, BOTH, RAISED, LEFT, CENTER, TOP, BOTTOM, Canvas, Scrollbar
from tkinter.ttk import Style, Entry
import json

def get_button(button,buttons):
    print(button)
    text = 'Updated'
    button['text'] = text
    print('Hello!')
    if button in buttons:
        print('Working')
        print(buttons.index(button))

class Get_Button_Pressed:
    def __init__(self, primary):
        #super().__init__(primary)
        self.primary = primary
        primary.geometry("200x200")
        primary.title("VME Digital to Analog Converter Control GUI")

        #Create a frame to hold directions and settings along with the save settings button.
        info_frame = Frame(primary, relief=RAISED, borderwidth=2)
        info_frame.pack(side='left', fill=BOTH, expand=True)

        buttons = []

        test_btn1 = Button(info_frame, text='Button 1', width=10, height=4, font='Helvetica 12 bold', bg = "red")
        test_btn1.configure(command=lambda button=test_btn1, buttons=buttons: get_button(button,buttons))
        buttons.append(test_btn1)
        test_btn1.pack(side="top")

        test_btn2 = Button(info_frame, text='Button 2', width=10, height=4, font='Helvetica 12 bold', bg = "green")
        test_btn2.configure(command=lambda button=test_btn2, buttons=buttons: get_button(button,buttons))
        #test_btn2 = Button(info_frame, text='Button 2', width=10, height=4, font='Helvetica 12 bold', command=self.greet, bg = "green")
        buttons.append(test_btn2)
        test_btn2.pack(side="top")

    #def greet(self,event):
    #    text = event.widget.cget('text')
    #    print(text)
    #    print('Pressed!')

    def greet(self):
        text = widget.cget('text')
        print(text)
        print('Pressed!')

    #def get_button():
    #    print('Hello!')

root = Tk()
my_gui = Get_Button_Pressed(root)
root.mainloop()
