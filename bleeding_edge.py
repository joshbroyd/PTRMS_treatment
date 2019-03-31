#At the least a list of the paths to each of the data files is needed.

from tkinter import (filedialog, Tk, Button, Label, StringVar, OptionMenu, 
    Entry, IntVar, Checkbutton)
import pandas as pd
import numpy as np
import sys, itertools, os
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def close():
    window.quit()
    window.destroy()

def getdatapaths():
    global paths
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
                tickbox.grid(row=5+n, column=m-1) #m+(len(channel_keys[0])/14 + 1)
                params[8][m].append(entry)
                tickboxes.append(tickbox)
    
    if params[2].get() == "Time series":
        for n in range(len(channels[0])):
            x = n//14 + 2
            y = n + 5 - (n//14 * 14)           
            entry = IntVar()
            tickbox = Checkbutton(window, text=channels[0][n], 
                        variable=entry)
            tickbox.grid(row=y, column=x, sticky="W")
            params[8][0].append(entry)
            tickboxes.append(tickbox)

    elif params[2].get() == "Mass scan":
        mz_channels = [n.split() for n in channels[0]]
        available_mz = ("m/z " + str(mz_channels[0][1]) + " to " 
            + str(mz_channels[-1][1]) + " are available")
        usr_inst = ("Input individual comma separated values\n"
            "or a range, hypen separated.")
        Label(window, text=available_mz).grid(row=5, column=2)
        Label(window, text=usr_inst).grid(row=7, column=2)
        params[8][0] = Entry(window)   
        params[8][0].insert("0", str(mz_channels[0][1]) + "-" 
            + str(mz_channels[-1][1]))
        params[8][0].grid(row=6, column=2)
    
def plot():
    if params[2].get() == "Time series":
        plot_time_series()
    elif params[2].get() == "Mass scan":
        plot_mass_scan()

def plot_time_series():
    print("plotting time series")

    _, ax1 = plt.subplots()
    if params[1].get() == "Absolute Time":
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))                    
                
    markercolours = itertools.cycle(['k','lightgreen','r','magenta',
                                     'midnightblue','darkorange'])
    linecolours = itertools.cycle(['grey','g','maroon','orchid', 'skyblue',
                                   'orange'])
    linestyles = itertools.cycle(['-', '--', ':', '-.'])
    
    print(params[1].get())
    
    all_channels = get_channels()

    ydata, ylabel, chosenchannels = get_ydata(params[1].get(), params[8][0], 
        all_channels[0])
    absolute_time = get_xdata("Absolute Time", params[3].get())[0]
    date = datetime.datetime.strftime(absolute_time[0], '%Y-%m-%d')
    
    reloffset = date + ' ' + params[3].get()
    
    xdata, xlabel = get_xdata(params[1], reloffset)   
    rel_time = get_xdata("Relative Time", 
        datetime.datetime.strftime(absolute_time[0], '%Y-%m-%d %H:%M:%S'))[0]
    
    
    cycles_perxmins = int(np.argmin([abs(element*60 - params[6].get())
                          for element in rel_time]))                          
    print("There are " + str(cycles_perxmins) + " cycles per " 
          + str(params[6].get()) + " seconds.")

    for index in range(len(ydata)):
        mc = next(markercolours)
        lc = next(linecolours)
        ls = next(linestyles)
        series_label = (chosenchannels[index] + ', ' + str(params[6].get())
        + ' second moving average')
        ysmooth = smooth(ydata[index], cycles_perxmins)
        ax1.plot(xdata, ydata[index], 'o',  ms=2, color=mc)
        ax1.plot(xdata, ysmooth, lw=2, color=lc, label=series_label,
                 linestyle=ls)
    
    plt.show()

def plot_mass_scan():
    print("plotting mass scan")
                
    yoption = params[0].get()
    all_channels = get_channels()
    channels = np.zeros(len(all_channels[0]), dtype=int)
    n_channels = convertparam(params[8][0].get())      
    n_channels = ["m/z " + str(i) for i in n_channels]
    for n in n_channels:
        channels[all_channels[0].index(n)] = 1
               
    y, ylabel, chosenchannels = get_ydata(yoption, channels, all_channels[0])
    print(ylabel, chosenchannels)
        
    files = paths
    lengths = []
    for f in files:
        data = np.asarray(pd.read_excel(f, sheet_name="Time   Cycle")
            ["Cycle number"])
        lengths.append(len(data))
 
    means = [[] for _ in range(len(files))]
    stddev = [[] for _ in range(len(files))]

    for m in range(len(chosenchannels)):
        count1 = 0
        for i in range(len(files)):
            count2 = count1 + lengths[i] 
            means[i].append(np.mean(y[m][count1:count2]))
            #Need to also calculate the standard error of the mean
            stddev[i].append(np.std(y[m][count1:count2]))
            count1 = count2
                
    ind = np.arange(len(chosenchannels))
    width = 2
        
    _, ax = plt.subplots()
        
    fccolours = itertools.cycle(['red', 'gray','green', 'blue', 'cyan',
        'magenta','lawngreen','darkorange','maroon','orchid'])
        
    for i in range(len(files)):
        files[i] = os.path.basename(files[i])[:-5]
        
    for i in range(len(files)):
        fc = next(fccolours)
        ax.bar(ind + ind*(len(files)*width) + width*i, means[i], width, 
            color=fc, alpha=.5, align='edge',
            error_kw=dict(ecolor='k', lw=1.5, capsize=3, capthick=1.5), 
            yerr=stddev[i], label=files[i])
        
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Mass/charge ratio")
    ax.set_xticks(ind + ind*(len(files)*width) + width*(len(files)/2))
    xlabels = [n[4:-2] for n in chosenchannels]
    ax.set_xticklabels(xlabels)
    ax.legend()
    plt.show()

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

def get_channels():
#Accepts a list of filenames and returns the full list of channel names (keys) 
#from the first excel file that's read in

    files = paths
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

def get_ydata(yoption, channels, channel_keys):
#Args are corresponding to raw counts or concentration and the list of channel 
#keys. Returns a list of lists of ydata from the corresponding channels and 
#yoption and the yaxis label. 
 
    chosenchannels = []
    for n in range(len(channels)):
        if channels[n] == 1:
            chosenchannels.append(channel_keys[n])
    
    y = [[] for _ in range(len(chosenchannels))]      
    for f in paths:
        data = pd.read_excel(f, sheet_name=yoption)
        for i in range(len(chosenchannels)):
            tmp = np.asarray(data[chosenchannels[i]])
            tmp = [0 if item == '#NV' else item for item in tmp]
            y[i].extend(tmp)
    
    if yoption == "Raw signal intensities":
            ylabel = "Raw signal intensity (cps)"
           
    elif yoption == "Concentration":
        ylabel = "Concentration (ppb)"
            
    elif yoption == "Instrument" or yoption == "Reaction conditions" :
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

    if reloffset != " " and xoption == "Relative Time":
        reloffset = datetime.datetime.strptime(reloffset.strip(),
                                               '%Y-%m-%d %H:%M:%S')
    
    if reloffset == " " and xoption == "Relative Time":
        sys.exit("Need to choose a starting time for relative time")
        
    x = []
    for f in paths:
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