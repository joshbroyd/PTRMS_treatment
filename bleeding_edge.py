#At the least a list of the paths to each of the data files is needed.

from tkinter import (filedialog, Tk, Button, Label, StringVar, OptionMenu, 
    Entry, IntVar, Checkbutton)
import pandas as pd
import numpy as np

def close():
        window.quit()
        window.destroy()

def getdatapaths():
    paths = list(filedialog.askopenfilenames(parent=window,
        title='Choose data files', 
        initialdir="/home/jgb509/Documents/CRM/Data"))
    params[5].insert("0", paths)
    
def getreadmepath():
    path = list(filedialog.askopenfilenames(parent=window,
        title='Choose readme file', 
        initialdir="/home/jgb509/Documents/CRM/Data"))
    params[4].insert("0", path)

def loaddata():
    global tickboxes
    tickboxes = []
    channels = get_channels()
    params[8] = [[] for _ in range(3)]

    for m in range(1, 3):
            for n in range(len(channels[m])):            
                entry = IntVar()
                tickbox = Checkbutton(window, text=channels[m][n],
                            variable=entry)
                tickbox.grid(row=5+n, column=m) #m+(len(channel_keys[0])/14 + 1)
                params[8][m].append(entry)
                tickboxes.append(tickbox)
    
    mz_channels = [n.split() for n in channels[0]]
    available_mz = ("m/z " + str(mz_channels[0][1]) + " to " 
            + str(mz_channels[-1][1]) + " are available")
    usr_inst = ("Input individual comma separated values\n"
            "or a range, hypen separated.")
    Label(window, text=available_mz).grid(row=5, column=1)
    Label(window, text=usr_inst).grid(row=7, column=1)
    params[8][0] = Entry(window)   
    params[8][0].insert("0", str(mz_channels[0][1]) + "-" 
            + str(mz_channels[-1][1]))
    params[8][0].grid(row=6, column=1)
    

def plot():
    print("hi")

def select_param():
    
    global params
    params = [[] for _ in range(9)]
    commands = [getdatapaths, getreadmepath, loaddata, plot, close]
    commandtext = ["Open files", "Readme file" , "Load data", "Plot", "Exit"] 
    for i in range(len(commandtext)):
        Button(window, text=commandtext[i], command=commands[i]).grid(row=i)

    Label(window, text="Moving average\nduration (seconds)").grid(row=3,
        column=2)
    Label(window, text="Absolute time (hh:mm) \nwhen relative time = 0").grid(
        row=4, column=2)
    
    params[7] = IntVar()
    Checkbutton(window, text="Calibration", 
        variable=params[7]).grid(column=2,row=2)

    paramsindex = range(3,7)
    rowindex = [3, 1, 0, 4]
    columnindex = [3, 1, 1, 3]
    for i in range(len(paramsindex)):
        params[paramsindex[i]] = Entry(window)
        params[paramsindex[i]].grid(column=columnindex[i], row=rowindex[i])

    menutext = ["y-axis format", "x-axis format",""]
    defaultoption = ["Concentration","Absolute time","Time series"]
    options = [["Concentration", "Raw signal intensities"], ["Cycle number", 
                "Absolute Time", "Relative Time"],["Time series", "Mass scan"]]
    for i in range(3):
        Label(window, text=menutext[i]).grid(row=i, column=2)
        params[i] = StringVar(window)
        params[i].set(defaultoption[i])
        OptionMenu(window, params[i], *options[i]).grid(row=i, column=3)
    
    print(params)

def get_channels():
#Accepts a list of filenames and returns the full list of channel names (keys) 
#from the first excel file that's read in

    files = params[5]
    all_channels = []
    sheetnames = ["Raw signal intensities", "Instrument",
                  "Reaction conditions"]                 
    for sheetname in sheetnames:
        data = pd.read_excel(files[0], sheet_name=sheetname)
        channels = [str(x) for x in list(data.keys())]
        all_channels.append(channels)
    
    #Need to check the channels of the other files in turn as well.

    return all_channels

# Converts the user input to parameter list
def convertparam(param):
    
    if '-' in param:
        i = param.index('-')
        MIN = float(param[:i])
        MAX = float(param[i+1:])
        STEP = 1.0
        NUM = int((MAX-MIN) / STEP+1)
        param = np.linspace(MIN, MAX, NUM)
        
    elif '-' not in param and ',' in param:
        param = [float(n) for n in param.split(',')]
    
    elif '-' not in param and ',' not in param:
        param = [float(param)]
       
    return param

def get_ydata(yoption, channels, channel_keys, sheetname):
#Args are corresponding to raw counts or concentration and the list of channel 
#keys. Returns a list of lists of ydata from the corresponding channels and 
#yoption and the yaxis label. 

    files = params[5]
    
    chosenchannels = []
    for n in range(len(channels)):
        if channels[n] == 1:
            chosenchannels.append(channel_keys[n])
    
    y = [[] for _ in range(len(chosenchannels))]      
    for f in files:
        data = pd.read_excel(f, sheet_name=sheetname)
        for i in range(len(chosenchannels)):
            tmp = np.asarray(data[chosenchannels[i]])
            tmp = [0 if item == '#NV' else item for item in tmp]
            y[i].extend(tmp)
    
    if sheetname == "Raw signal intensities" or sheetname == "Concentration":
        if yoption == "Raw signal intensities":
            ylabel = "Raw signal intensity (cps)"
           
        elif yoption == "Concentration":
            ylabel = "Concentration (ppb)"
            
    elif sheetname == "Instrument" or sheetname == "Reaction conditions" :
        ylabel = "ARB"
    
   # print(sheetname, yoption)    
   # print(ylabel)     
    return y, ylabel, chosenchannels

def smooth(y, box_pts):
#Smooth function to return the moving average of length box_pts from list y.     
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def get_xdata(xoption, reloffset):
#Returns xdata list and xlabel depending on xoption.

    files = params[5]
    if reloffset != " " and xoption == "Relative Time":
        reloffset = datetime.datetime.strptime(reloffset.strip(),
                                               '%Y-%m-%d %H:%M:%S')
    
    if reloffset == " " and xoption == "Relative Time":
        sys.exit("Need to choose a starting time for relative time")
        
    x = []
    for f in files:
        data = pd.read_excel(f, sheet_name="Time   Cycle")
        if xoption == "Cycle number":
            data = np.asarray(data[xoption])
            offset = len(x)
            x.extend([item + offset for item in data])
            xlabel = xoption + " (ARB)"
                    
        elif xoption == "Absolute Time":
            tmp = []
            for time in range(len(data[xoption])):
                tmp.append(data[xoption][time].to_pydatetime())
            x.extend(tmp)
            xlabel = xoption + " (hh:mm)"
                    
        elif xoption == "Relative Time":          
            tmp = []
            for time in range(len(data["Absolute Time"])):
                tmp.append(data["Absolute Time"][time].to_pydatetime())
            
            tmp = [(item-reloffset).total_seconds()/60.0 for item in tmp]
            x.extend(tmp)
            xlabel = xoption + " (mins)"
                   
    return x, xlabel

window = Tk()
select_param()
window.mainloop()