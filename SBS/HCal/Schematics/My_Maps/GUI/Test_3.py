from tkinter import *

root = Tk(); root.geometry( "225x300+530+180" )
buttons = {}

def makeChoice( event ):
    global buttons
    print(event,'  ',event.widget,'  ',buttons[ event.widget ])
    
def createBoard():
    global buttons
    buttonNum = 0
    x = -225
    y = 65
    for b in range( 9 ):
        for b2 in range( 9 ):
            button = Button( root, text = " ", font = "Courier 9", width = 2, bd = 3); button.place( relx = 1, x = x, y = y )
            buttons[ button ] = buttonNum
            buttonNum += 1
            button.bind( "<Button-1>", makeChoice )
            x += 25
        x = -225
        y += 26

createBoard()
mainloop()
