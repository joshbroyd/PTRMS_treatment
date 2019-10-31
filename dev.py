from tkinter import (Tk, filedialog, Entry, Button, IntVar,
    Checkbutton, OptionMenu, StringVar, Label)
import pandas as pd
import numpy as np
import os

def smooth(y, box_pts):
###Returns the moving average of length box_pts from list y.     
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def GetExcelFiles():
###Creates the filedialog for the excel files
    initialdirpath = "/Measurements/PTR-MS/Data"
    f = filedialog.askopenfilenames(initialdir=initialdirpath)
    ExcelFiles.delete("0", "-1")
    ExcelFiles.insert("0", f)    

def GetReadmeFile():
###Creates the filedialog for the readme
    initialdirpath = "/Measurements/PTR-MS/Data"
    openfilename = filedialog.askopenfilename(initialdir=initialdirpath)
    ReadmeFile.delete("0", "-1")
    ReadmeFile.insert("0", openfilename)


def ConvertParam(param):
###Converts user input to list    
    if '-' in param and ',' not in param:
        i = param.index('-')
        MIN = int(param[:i])
        MAX = int(param[i+1:])
        param = range(MIN, MAX + 1)
    elif '-' in param and ',' in param:
        param = param.split(',')
        for m in range(len(param)):
            if '-' in param[m]:
                i = param[m].index('-')
                MIN = int(param[m][:i])
                MAX = int(param[m][i+1:])
                param[m] = range(MIN, MAX + 1)
            else:
                param[m] = int(param[m])
    elif '-' not in param and ',' in param:
        param = [int(n) for n in param.split(',')]
    elif '-' not in param and ',' not in param:
        param = [int(param)]
    return param

def GetKeys():
###Returns list of Sheetnames and Keys
    Sheetnames = ['Time   Cycle', 'Raw signal intensities', 'Concentration',
                  'Instrument', 'Reaction conditions']
    Keys = [[] for _ in range(len(Sheetnames))]
    for f in ExcelFiles.get().split():
        for i in range(len(Keys)):
            Keys[i].extend(list(pd.read_excel(f, Sheetnames[i]).keys()))
    for i in range(len(Keys)):
        Keys[i] = list(set(Keys[i]))
        Keys[i].sort()
    Keys = Keys[:2] + Keys[3:]
    Keys[2] = [Keys[2][0]] + Keys[2][3:9] + Keys[2][10:]
    Keys[3] = Keys[3][6:]
    Keys = [Keys[0]] + [Keys[2] + Keys[3]] + [Keys[1]]
    return Sheetnames, Keys

def LoadData():
    Sheetnames, Keys = GetKeys()
    if Tickboxes != []:
        for t in Tickboxes:
            t.destroy()
    if GUILabels != []:
        for l in GUILabels:
            l.destroy()
        ChannelEnts[-1].destroy()

    for i in range(len(Keys[1])):
        x = i//14
        y = i + 5 - (i//14 * 14)
        Ent = IntVar()
        Tickbox = Checkbutton(root, text=Keys[1][i], variable=Ent)
        Tickbox.grid(column=x, row=y, sticky="NESW")
        ChannelEnts.append(Ent)
        Tickboxes.append(Tickbox)
    
    Label(root, text="Quadrupole channel\nformat").grid(column=2, row=0, sticky="NSEW")
    QuadFormat = StringVar()
    QuadFormat.set(Sheetnames[2])
    OptionMenu(root, QuadFormat, *Sheetnames[1:3]).grid(column=3, row=0, sticky="NSEW")
    
    if PlotFormat.get() == "Time series":
        for i in range(len(Keys[2])):
            x = i//14 + 1
            y = i + 5 - (i//14 * 14)
            Ent = IntVar()
            Tickbox = Checkbutton(root, text=Keys[2][i], variable=Ent)
            Tickbox.grid(column=x, row=y, sticky="NESW")
            ChannelEnts.append(Ent)
            Tickboxes.append(Tickbox)

        TimeLabels = ["X axis format", "Moving average\nduration (seconds)"]
        for i in range(len(TimeLabels)):
            l = Label(root, text=TimeLabels[i])
            Tickboxes.append(l)
            l.grid(column=2, row=i+1, sticky="NSEW")

        XaxisFormat = StringVar()
        XaxisFormat.set(Keys[0][1])
        XaxisMenu = OptionMenu(root, XaxisFormat, *Keys[0])
        XaxisMenu.grid(column=3, row=1, sticky="NSEW")
        Tickboxes.append(XaxisMenu)
        
        MovingAvEnt = Entry()
        MovingAvEnt.insert("0", "120")
        MovingAvEnt.grid(column=3, row=2, sticky="NSEW")
        Tickboxes.append(MovingAvEnt)

    
    elif PlotFormat.get() == "Mass scan":
        m = [Keys[2][0].split()[1], Keys[2][-1].split()[1]]
        usr_inst = ["m/z " + m[0] + " to " + m[1] + " are available",
                    ("Input individual comma separated values\n"
                     "or a range, hypen separated.")]
        MassScanEnt = Entry()
        ChannelEnts.append(MassScanEnt)   
        MassScanEnt.insert("0", m[0] + "-" + m[1])
        MassScanEnt.grid(row=6, column=1)
        for i in range(len(usr_inst)):
            l = Label(root, text=usr_inst[i])
            GUILabels.append(l)
            l.grid(row=(i+1)*2 + 3, column=1, sticky="NESW")

def PlotData():
    print("TO DO")

root = Tk()

Tickboxes = []
GUILabels = []
ChannelEnts = []

Button(root, text="Choose excel files", command=GetExcelFiles).grid(column=0, row=0, sticky="NSEW")
ExcelFiles = Entry(root)
ExcelFiles.grid(column=1, row=0, sticky="NSEW")
ExcelFiles.insert("0", "")

Button(root, text="Choose Readme", command=GetReadmeFile).grid(column=0, row=1, sticky="NSEW")
ReadmeFile = Entry(root)
ReadmeFile.grid(column=1, row=1, sticky="NSEW")
ReadmeFile.insert("0", "")

Button(root, text="Load Data", command=LoadData).grid(column=0, row=2, sticky="NSEW")

PlotFormat = StringVar()
PlotFormat.set("Time series")
OptionMenu(root, PlotFormat, *["Time series", "Mass scan"]).grid(column=1, row=2,sticky="NESW")

Button(root, text="Plot", command=PlotData).grid(column=0, row=3, sticky="NSEW")
Button(root, text="Exit", command=root.quit).grid(column=0, row=4, sticky="NSEW")
root.mainloop()
root.destroy()