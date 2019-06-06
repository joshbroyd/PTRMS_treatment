from tkinter import (Tk, Frame, Button, filedialog, Entry, IntVar, Checkbutton)
import numpy as np
import matplotlib.pyplot as plt

class Application(Frame):

    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.grid()
        self.createWidgets()
        

    def createWidgets(self):
        
        ctext = ["Add Excel Files", "Add Project File", 
                 "Add Broadband Spectra Folder",
                 "Add OHN2 Spectra Folder"] 
        dir1 = "/home/jgb509/Documents/CRM/Data"
        dir2 = "/home/jgb509/Documents/CRM/Spectroscopy"
        idirs = [dir1, dir1, dir2, dir2]

        commands = [getpath(idirs[n]) for n in range(len(ctext))]
        self.paths = [[] for _ in range(len(ctext))]

        for i in range(len(ctext)):
            Button(self, text=ctext[i], 
                   command=commands(i)).grid(row=i, sticky="NSEW")
            self.paths[i] = Entry(self)
            self.paths[i].grid(row=i, column=1, sticky="NSEW")
        
        Button(self, text="Exit", command=self.quit).grid(sticky="NSEW")

#PROBLEM: Cannot have a list of inital dirs without the filedialog opening.

def getpath(i):
    filepath = filedialog.askopenfilename()#initialdir=app.idirs[i])

  #  app.paths[i].insert("0", filepath) 



root = Tk()
app = Application(master=root)
app.mainloop()
root.destroy()
