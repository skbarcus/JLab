from tkinter import Tk, RIGHT, BOTH, RAISED, LEFT, CENTER
from tkinter.ttk import Frame, Button, Style

class Example(Frame):

    def __init__(self):
        super().__init__()

        self.initUI()


    def initUI(self):

        self.master.title("HCal")
        self.style = Style()
        self.style.theme_use("default")

        frame1 = Frame(self, relief=RAISED, borderwidth=2)
        frame1.pack(fill=BOTH, expand=True)

        frame2 = Frame(self, relief=RAISED, borderwidth=2)
        frame2.pack(fill=BOTH, expand=True)

        self.pack(fill=BOTH, expand=True)

        testButton1 = Button(frame1, text="Test 1")
        testButton1.pack(side="left", padx=5, pady=5)

        testButton2 = Button(frame2, text="Test 2")
        testButton2.pack(side="left", padx=5, pady=5)

        testButton3 = Button(self, text="Test 3")
        testButton3.pack(side="left", padx=5, pady=5)

        testButton4 = Button(self, text="Test 4")
        testButton4.pack(side="left", padx=5, pady=5)

        testButton5 = Button(self, text="Test 5")
        testButton5.pack(side="left", padx=5, pady=5)

        testButton6 = Button(self, text="Test 6")
        testButton6.pack(side="left", padx=5, pady=5)

        closeButton = Button(self, text="Close")
        closeButton.pack(side="left", padx=5, pady=5)
        okButton = Button(self, text="OK")
        okButton.pack(side="left")


def main():

    root = Tk()
    root.geometry("600x800")
    #root.geometry("600x800+300+300")
    app = Example()
    root.mainloop()


if __name__ == '__main__':
    main()
