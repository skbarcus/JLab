class MyApp(object):
    def __init__(self):
        self.root = Tk()
        self.root.wm_title("ImagePro")

        #Original
        original = Image.open(infile)
        (w, h) = (original.size[0], original.size[1])
        tkpi = ImageTk.PhotoImage(original)
        label = Label(self.root, image=tkpi)
        label.grid(row =0, column=0, padx=5,pady=5)

        img = original.copy().convert("L")
        tkpi2 = ImageTk.PhotoImage(img)
        label = Label(self.root, image=tkpi2)
        label.grid(row =0, column=1, padx=5,pady=5)

        Label(self.root, text = "Original").grid(row=1, column=0)
        Label(self.root, text = "Modified").grid(row=1, column=1)

        self.buttonframe = Frame(self.root)
        self.buttonframe.grid(row=2, column=0, columnspan=2)        
        Button(self.buttonframe, text = "Brighten").grid(row=0, column=0)
        Button(self.buttonframe, text = "Darken").grid(row=0, column=1)
        Button(self.buttonframe, text = "Warm").grid(row=0, column=2)
        Button(self.buttonframe, text = "Cool").grid(row=0, column=3)

        self.root.mainloop()
