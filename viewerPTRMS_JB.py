#!/usr/bin python3

from tkinter import (Tk, Frame, Button, filedialog, Entry, IntVar, Checkbutton,
Label, StringVar, OptionMenu)
#mpl.use('qt5agg')
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import sys, itertools, datetime, bisect, string
from scipy.stats import linregress
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse.linalg import spsolve
import scipy.sparse as sparse
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import pandas as pd

#Changes the font and fontsize of the graphs
import os
import matplotlib as mpl
dir_path = os.path.dirname(os.path.realpath(__file__))
if __name__ == "__main__":

    fontsize = 25
    params = {'backend':'QT5Agg',
            'text.latex.preamble':['\\usepackage{gensymb}'],
            'text.usetex': True,
            'axes.labelsize':fontsize,
            'axes.titlesize':fontsize,
            'font.size':fontsize,
            'legend.fontsize':fontsize-5,
            'xtick.labelsize':fontsize,
            'ytick.labelsize':fontsize,
            'xtick.direction':'in',
            'xtick.major.size':10,
            'xtick.major.width':2,
            'xtick.minor.size':5,
            'xtick.minor.width':2,
            'ytick.direction': 'in',
            'ytick.major.size':10,
            'ytick.major.width':2,
            'ytick.minor.size':5,
            'ytick.minor.width':2,
            'lines.markersize':10,
            'lines.linewidth':3,
            'axes.linewidth':2,
            'axes.grid':True,
            'figure.figsize': [16, 9],
            'font.family':'serif',
            'font.serif':'Times'}
    mpl.rcParams.update(params)

#plt.rc('text', usetex=True)

class Application(Frame):

    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.grid()
        self.createWidgets()
        
    def createWidgets(self):
        r, c = [0, 1, 0, 2], [0, 0, 4, 4]
        ctext = ["Add Excel Files", "Add Readme", 
                 "Add Broadband Spectra Folder",
                 "Add OHN2 Spectra Folder"]
        commands = [self.getdatapath, self.getreadmepath, 
                    self.getbroadpath, self.getOHN2path]
        self.paths = [[] for _ in range(len(ctext))]
        #self.paths[0] = Excel data path list
        # "[1] = Readme path
        # "[2] = Broadband spectra folder
        # "[3] = OHN2 spectra folder            
        for i in range(len(ctext)):
            Button(self, text=ctext[i], command=commands[i]).grid(
                row=r[i], column=c[i], sticky="NSEW")
            self.paths[i] = Entry(self)
            self.paths[i].grid(row=r[i], column=c[i]+1)
        
        self.params = [[] for _ in range(9)]
        #self.params[0] = Absolute time when relative time = 0
        # "[1] = Moving average duration
        # "[2] = Broadband lines
        # "[3] = OHN2 lines
        # "[4] = ylimits
        # "[5] = y-axis format
        # "[6] = x-axis format
        # "[7] = time series or mass scan
        # "[8] = calibration
        #baselineparams
        # "[0] = lambda parameter for baseline correction
        # "[1] = p parameter for baseline correction
        # plot and use baseline correction
        self.baselineparams = [[] for _ in range(2)]
        labels = ["Time to take\nbackground", "Calibration factor(s) (,sep)"]
        defaults = ["hh:mm:ss - hh:mm:ss", 1]
        r, c = [0, 1], [7, 7]
        for i in range(len(labels)):
            Label(self, text=labels[i]).grid(column=c[i]-1, row=r[i])
            self.baselineparams[i] = Entry(self)
            self.baselineparams[i].grid(column=c[i],row=r[i])
            self.baselineparams[i].insert("0", defaults[i])

        Label(self, text="Starting annotation height").grid(column=6, row=3, sticky="NSEW")
        self.arrow_heightent = Entry(self)
        self.arrow_heightent.grid(column=7, row=3, sticky="NSEW")
        self.arrow_heightent.insert("0", "100")

        Label(self, text="Annotation height offset").grid(column=6, row=4, sticky="NSEW")
        self.arrow_heightoffsetent = Entry(self)
        self.arrow_heightoffsetent.grid(column=7, row=4, sticky="NSEW")
        self.arrow_heightoffsetent.insert("0", "3")

        r, c = [3, 4, 1, 3, 4], [3, 3, 5, 5, 5]
        labeltext = ["Absolute time (hh:mm:ss) \nwhen trel=0",
                     "Moving average\nduration (seconds)",
                     "Broadband lines \n (O peak at 777.25nm and 844.66nm):",
                     ("OHN2 lines (OH peak at 308.92nm, \nN2 peaks at " 
                      "336.30nm, 357.56nm):"), "y limits (l-axis, r-axis,\nspec-axis, hypen delimited)"]
        defaults = ["hh:mm:ss", 120.0, "777.25, 844.66", "308.92, 336.30, 357.56", 
                    "ymin - ymax"]
        for i in range(len(labeltext)):
            Label(self, text=labeltext[i]).grid(column=c[i]-1, row=r[i])
            self.params[i] = Entry(self)
            self.params[i].grid(column=c[i],row=r[i])
            self.params[i].insert("0", defaults[i])
        
        menutext = ["Quadrupole channel \n format", "x-axis format",""]
        defaultoption = ["Concentration", "Absolute Time", "Time series"]
        options = [["Concentration", "Raw signal intensities", "Normalised signal intensity"], 
                   ["Cycle number", "Absolute Time", "Relative Time"],
                   ["Time series", "Mass scan"]]
        for i in range(len(menutext)):
            Label(self, text=menutext[i]).grid(row=i, column=2, sticky="NESW")
            self.params[i+5] = StringVar(self)
            self.params[i+5].set(defaultoption[i])
            OptionMenu(self, self.params[i+5], *options[i]).grid(
                row=i, column=3, sticky="NESW")
        self.params[8] = IntVar()
        Checkbutton(self, text="Calibration", variable=self.params[8]).grid(
            column=2, row=2, sticky="NSEW")
        self.usebaseline = IntVar()
        Checkbutton(self, text="Use baseline correction\n& calibration factors", variable=self.usebaseline).grid(
            column=7, row=2, sticky="NSEW")
        self.logy = IntVar()
        Checkbutton(self, text="Logy", variable=self.logy).grid(
            column=6, row=2, sticky="NSEW")
        self.EN = IntVar()
        Checkbutton(self, text="Plot E/N?", variable=self.EN).grid(
            column=1, row=2, sticky="NSEW")
        self.yaxis2 = IntVar()
        Checkbutton(self, text="Use 2nd yaxis?", variable=self.yaxis2).grid(
            column=1, row=3, sticky="NSEW")
        self.specden = IntVar()
        Checkbutton(self, text="Plot species density?", variable=self.specden).grid(
            column=1, row=4, sticky="NSEW")
        Button(self, text="Load data", command=self.loaddata).grid(
            row=2, column=0, sticky="NSEW")
        Button(self, text="Plot", command=plot).grid(
            row=3,column=0, sticky="NSEW")
        Button(self, text="Exit", command=self.quit).grid(
            row=4, column=0, sticky="NSEW")

    def getdatapath(self):
        path = list(filedialog.askopenfilenames(
            title="Choose data files",
            initialdir="/home/jgb509/Documents/Measurements/PTR-MS/Data"))
        mpl.rcParams["savefig.directory"] = os.path.dirname(os.path.abspath(path[0]))#
        mpl.rcParams["savefig.format"] = "eps"
        self.paths[0].insert("0", path)

    def getreadmepath(self):
        path = filedialog.askopenfilename(
            title="Choose readme file",
            initialdir="/home/jgb509/Documents/Measurements/PTR-MS/Data")
        self.paths[1].insert("0", path)

    def getbroadpath(self):
        path = filedialog.askdirectory(
            title='Choose broadband spectroscopy folder', 
            initialdir="/home/jgb509/Documents/Measurements/Spectroscopy")
        self.paths[2].insert("0", path)    

    def getOHN2path(self):
        path = filedialog.askdirectory(
            title='Choose OHN2 spectroscopy folder', 
            initialdir="/home/jgb509/Documents/Measurements/Spectroscopy")
        self.paths[3].insert("0", path)

    def loaddata(self):
        all_channels, transmission = get_channels()
        self.tickboxes = []
        self.channels = []
        # first 23 entries are instrumental conditions.
        
        for m in range(23):
            x = m//14
            y = m + 5 - (m//14 * 14)
            entry = IntVar()
            tickbox = Checkbutton(self, text=all_channels[m], 
                                    variable=entry)
            tickbox.grid(row=y, column=x, sticky="NESW")
            self.channels.append(entry)
            self.tickboxes.append(tickbox)
        
        if self.params[7].get() == "Time series":
            for n in range(23, 23+len(all_channels[23:])):
                x = (n-23)//14 + 2
                y = (n-23) + 5 - ((n-23)//14 * 14)           
                entry = IntVar()
                tickbox = Checkbutton(self, text=all_channels[n], 
                                      variable=entry)
                tickbox.grid(row=y, column=x, sticky="NESW")
                self.channels.append(entry)
                self.tickboxes.append(tickbox)
        
        elif self.params[7].get() == "Mass scan":
            mz_channels = [n.split() for n in all_channels[23:]]
            available_mz = ("m/z " + str(mz_channels[0][1]) + " to " 
                + str(mz_channels[-1][1]) + " are available")
            usr_inst = ("Input individual comma separated values\n"
                        "or a range, hypen separated.")
            Label(self, text=available_mz).grid(row=5, column=2, sticky="NESW")
            Label(self, text=usr_inst).grid(row=7, column=2, sticky="NESW")
            self.channels.append(Entry(self))   
            self.channels[23].insert("0", str(mz_channels[0][1]) + "-" 
                + str(mz_channels[-1][1]))
            self.channels[23].grid(row=6, column=2)

def get_channels():
#Accepts a list of filenames and returns the full list of channel names (keys) 
#from the first excel file that's read in

    all_channels = []
    sheetnames = ["Instrument","Reaction conditions"]
    paths = app.paths[0].get().split()
    for sheetname in sheetnames:
        data = pd.read_excel(paths[0], sheet_name=sheetname)
        channels = [str(x) for x in list(data.keys())]
        all_channels.extend(channels)
        #use set to have a unique list of channel names. Need to check each
        #individual file to check whether there is a channel of that name.
    
    tmp = []
    transmission = []
    for path in paths:
        data1 = pd.read_excel(path, sheet_name="Raw signal intensities")
        channels = [str(x) for x in list(data1.keys())]
        data = pd.read_excel(paths[0], sheet_name="Untitled (root)")
        tmp1 = np.asarray(data[[str(x) for x in list(data.keys())][17]])
        for c in range(len(channels)):
#time series channels are "m/z 21.00 ch0", mass scan channels are "m/z 21.0"
#This changes the internal keynames to ints, i.e. "m/z 21" they are changed
#depending on mass scan or time series later.
            tmp2 = channels[c].split()
            tmp3 = tmp2[0] + " " + str(int(float(tmp2[1])))
            tmp4 = tmp3 + " " + str(tmp1[18+c])
            tmp.append(tmp3)
            transmission.append(tmp4)
    transmission = list(set(transmission))
    transmission.sort()
    tmp = list(set(tmp))
    tmp.sort()
    all_channels.extend(tmp)
    return all_channels, transmission

# Converts the user input to parameter list
def convertparam(param):
    
    print(param)
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

def get_xdata(xoption, filenames):
#Returns xdata list, xlabel depending on xoption, date and cycles per xmins

    cycle = [[] for _ in range(len(filenames))]
    absolute = [[] for _ in range(len(filenames))]
    rel_time = [[] for _ in range(len(filenames))]
    tmp2 = [[] for _ in range(len(filenames))]
    cycles_perxmins = []

    for i in range(len(filenames)):
        for f in filenames[i]:
            data = pd.read_excel(f, sheet_name="Time   Cycle")

            cdata = np.asarray(data["Cycle number"])
            offset = len(cdata)
            cycle.extend([item + offset for item in cdata])

            tmp = []
            for time in range(len(data["Absolute Time"])):
                tmp.append(data["Absolute Time"][time].to_pydatetime())
            absolute[i].extend(tmp)
                
    for i in range(len(filenames)):
        tmp = []
        
#This first section just calculates cycles_per_xmins
        date = datetime.datetime.strftime(absolute[i][0], '%Y-%m-%d')
        print(date)
        f = filenames[i][0]
        data = pd.read_excel(f, sheet_name="Time   Cycle")
        reloffset = absolute[i][0]
        for time in range(len(absolute[i])):
            tmp.append(absolute[i][time])
        tmp = [(item-reloffset).total_seconds()/60.0 for item in tmp]
        tmp2[i].extend(tmp)
##Have to calculate cycles per xmins before using the reloffset otherwise cocks up ysmooth.
        cycles_perxmins.append(int(np.argmin([abs(element*60 - float(app.params[1].get())) for element in tmp2[i]])))
        tmp = []

        if app.params[0].get() != "hh:mm:ss":
            date = datetime.datetime.strftime(absolute[i][0], '%Y-%m-%d')
            reloffset = date + ' ' + app.params[0].get()
            print("Relative time starts at " + reloffset)
            reloffset = datetime.datetime.strptime(reloffset.strip(),
                                               '%Y-%m-%d %H:%M:%S')
        elif app.params[0].get() == "hh:mm:ss":
            date = datetime.datetime.strftime(absolute[i][0], '%Y-%m-%d')
            f = filenames[i][0]
            data = pd.read_excel(f, sheet_name="Time   Cycle")
            reloffset = absolute[i][0]   
            print("Relative time starts at " + datetime.datetime.strftime(reloffset,'%H:%M:%S'))
        for time in range(len(absolute[i])):
            tmp.append(absolute[i][time])
        tmp = [(item-reloffset).total_seconds()/60.0 for item in tmp]
        rel_time[i].extend(tmp)
                
    if xoption == "Cycle number":
        xlabel = xoption + " (ARB)"
        x = cycle
    elif xoption == "Absolute Time":
        xlabel = xoption + " (hh:mm)"
        x = absolute
    elif xoption == "Relative Time":
        xlabel = xoption + " (mins)"
        x = rel_time

    return x, absolute, xlabel, cycles_perxmins, date

def get_ydata(allchannels, channel_keys, transmission, filenames):
#Args are the list of channels & channel keys. Returns a list of lists of 
# ydata from the corresponding channels and the labels. 

#Need to know which sheetname to use for which channels. Needs to be adaptive
#so it will still plot if a mixture of sheetnames is chosen.

#Need to have a list of channel_keys that the user wants to plot. Needs to be organised 
#so that if just quadrupole channels want to be plotted the list of lists should look 
#like: [[],[],[usrchosen1, usrchosen2]] etc

    ydata = []
    rawsignal = []
    Udrift, Tdrift, pdrift = [], [], []
    ENvalues = [Udrift, Tdrift, pdrift]
    linelabels = []
    chosenkeychannels = [[] for _ in range(3)]

    InstrumentChannels = channel_keys[:14]
    ReacCondChannels = channel_keys[14:23]
    QuadChannels = channel_keys[23:]
    keychannels = [InstrumentChannels, ReacCondChannels, QuadChannels]

    InstrumentChannels = allchannels[:14]
    ReacCondChannels = allchannels[14:23]
    QuadChannels = allchannels[23:]
    chosenchannels = [InstrumentChannels, ReacCondChannels, QuadChannels]

    if app.params[5].get() == "Raw signal intensities":
        ylabel = "Raw signal intensity (cps)"
        sheetnames = ["Instrument","Reaction conditions", app.params[5].get()]
           
    elif app.params[5].get() == "Concentration":
        ylabel = "Concentration (ppbv)"
        sheetnames = ["Instrument","Reaction conditions", app.params[5].get()]
            
    elif app.params[5].get() == "Instrument" or app.params[5].get() == "Reaction conditions" :
        ylabel = "ARB"
        sheetnames = ["Instrument","Reaction conditions", app.params[5].get()]
    
    elif app.params[5].get() == "Normalised signal intensity":
        if app.usebaseline.get() == 1:
            ylabel = "Concentration (ppbv)"
        else:
            ylabel = "Normalised signal intensity (ncps)"
        sheetnames = ["Instrument","Reaction conditions", "Raw signal intensities"]
    
    for m in range(3):
        for n in range(len(chosenchannels[m])):
            if chosenchannels[m][n] == 1:
                chosenkeychannels[m].append(keychannels[m][n])
                linelabels.append(keychannels[m][n])
    
    chosenkeychannels[2] = list(reversed(chosenkeychannels[2]))
    linelabels = list(reversed(linelabels))
    print(chosenkeychannels[2])

#Separate into instrumental/reaction condition keys and raw sig/conc keys
#This part only looks at the instrument and reaction conditions sheets
    for n in range(len(sheetnames[:2])):
        for i in range(len(chosenkeychannels[n])):
            tmp2 = []
            for f in filenames[i]:
                data = pd.read_excel(f, sheet_name=sheetnames[n])
                tmp1 = np.asarray(data[chosenkeychannels[n][i]])
                tmp1 = [0 if item == '#NV' else item for item in tmp1]
                tmp2.extend(tmp1)
            ydata.append(tmp2)

    filenames = [filenames for _ in range(len(linelabels))]
    files_to_add = [[] for _ in range(len(linelabels))]
#This part only looks at the concentration and raw signal intensities sheets.
    for i in range(len(chosenkeychannels[2])):
        tmp2 = []
        tmp4 = []
        for f in filenames[i]:
#Need to know whether this is a mass scan or time series for the chosenkeychannel names.
#Time series: "m/z 21.00 ch0", mass scan: "m/z 21.0"
            if pd.read_excel(f, sheet_name="Untitled (root)")["Title"][0] == "SCAN data ":
                tmp3 = chosenkeychannels[2][i] + ".0"
            elif pd.read_excel(f, sheet_name="Untitled (root)")["Title"][0] == "MID data ":
                tmp3 = chosenkeychannels[2][i] + ".00"
            data = pd.read_excel(f, sheet_name=sheetnames[2])
            rawsigdata = pd.read_excel(f, sheet_name="Raw signal intensities")
            tmpkeys = [str(x) for x in list(data.keys())]
            if tmp3.split()[1] in [z.split()[1] for z in tmpkeys]:
                files_to_add[i].append(f)
#Fails when two different channels are selected with one not existing in one set of files.
#If the first file it selects has more datapoints then when it tries ot plot the second
#channel it fails. Need to have a list of filenames for each channel for the x-axis to
#pull out data from.
#Need to also have n sets of x-axis data lists.
                for key1 in tmpkeys:
                    if tmp3.split()[1] == key1.split()[1]:
                        tmp1 = np.asarray(data[key1])
                        tmp1 = [0 if item == '#NV' else item for item in tmp1]
                        tmp2.extend(tmp1)
                        tmp1 = np.asarray(rawsigdata[key1])
                        tmp1 = [0 if item == '#NV' else item for item in tmp1]
                        tmp4.extend(tmp1)
            ENdata = pd.read_excel(f, sheet_name="Reaction conditions")
            Udrift.extend([0 if item == "#NV" else item for item in np.asarray(ENdata["Udrift"])])
            pdrift.extend([0 if item == "#NV" else item for item in np.asarray(ENdata["pdrift"])])
            Tdrift.extend([0 if item == "#NV" else item+273.15 for item in np.asarray(ENdata["Tdrift"])])

        ydata.append(tmp2)
        rawsignal.append(tmp4)
#IMPORTANT: If there is not one of this chosen key in the file, the corresponding x-values
#need to be removed. Would be very easy to remove the file from the GUI list.

#This section below can screw around with the order of the files.
    for i in range(len(files_to_add)):
        files_to_add[i] = list(set(files_to_add[i]))
        filenames[i] = files_to_add[i]
    
    for i in range(len(filenames)):
        filenames[i] = filenames_sort(filenames[i])

    if app.params[5].get() == "Normalised signal intensity":
        for i in range(len(chosenkeychannels[2])):
            H3Odata = []
            for f in filenames[i]:
#Need to know whether this is a mass scan or time series for the chosenkeychannel names.
#Time series: "m/z 21.00 ch0", mass scan: "m/z 21.0"
                if pd.read_excel(f, sheet_name="Untitled (root)")["Title"][0] == "SCAN data ":
                    tmp3 = chosenkeychannels[2][i] + ".0"
                elif pd.read_excel(f, sheet_name="Untitled (root)")["Title"][0] == "MID data ":
                    tmp3 = chosenkeychannels[2][i] + ".00"
                data = pd.read_excel(f, sheet_name=sheetnames[2])
                tmpkeys = [str(x) for x in list(data.keys())]
                for key1 in tmpkeys:
                    if int(float(key1.split()[1])) == 21:
                        tmp1 = np.asarray(data[key1])
                        tmp1 = [0 if item == '#NV' else item for item in tmp1]
                        H3Odata.extend(tmp1)
        
            for n in range(len(linelabels)):
                if int(linelabels[n].split()[1]) == 21:
                    tmp = (1e6*0.7)/(np.array(ydata[i])*500)
                for m in range(len(transmission)):
                    if int(linelabels[n].split()[1]) == int(transmission[m].split()[1]) and int(linelabels[n].split()[1]) != 21:
                        tmp = (1e6*np.array(ydata[i])*0.7)/(np.array(H3Odata)*500*float(transmission[m].split()[2]))
            ydata[i] = tmp

    print("Chosen channels:", linelabels)
    return ydata, ylabel, linelabels, rawsignal, filenames, ENvalues

def smooth(y, box_pts):
#Smooth function to return the moving average of length box_pts from list y.     
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='valid')
    return y_smooth

def plot():
    if app.params[7].get() == "Time series":
        plot_time_series()
      #  for tickbox in app.tickboxes:
      #      tickbox.destroy()

    elif app.params[7].get() == "Mass scan":
        plot_mass_scan()
     #   for tickbox in app.tickboxes:
     #       tickbox.destroy()

def plot_mass_scan():
    print("plotting mass scan")
        
#Is it worth having a readme just for mass scans that indicates what times
#to analyse?

    all_channels, transmission = get_channels()
    print(all_channels)
    channels = np.zeros(len(all_channels), dtype=int)
    n_channels = convertparam(app.channels[23].get())

    n_channels = ["m/z " + str(i) for i in n_channels]
    print(n_channels)
    for n in n_channels:
        channels[all_channels.index(n)] = 1
               
    y, ylabel, chosenchannels, rawsignal, filenames = get_ydata(channels, all_channels, transmission)
    print(ylabel, chosenchannels)
    lengths = []

    paths = app.paths[0].get().split()
    paths.sort()
    print(paths)
    for f in paths:
        data = np.asarray(pd.read_excel(f, sheet_name="Time   Cycle")
            ["Cycle number"])
        lengths.append(len(data))
 
    means = [[] for _ in range(len(paths))]
    stddev = [[] for _ in range(len(paths))]
    stderr = [[] for _ in range(len(paths))]

    for m in range(len(chosenchannels)):
        count1 = 0
        for i in range(len(paths)):
            count2 = count1 + lengths[i] 
            means[i].append(np.mean(y[m][count1:count2]))
            stddev[i].append(np.std(y[m][count1:count2]))
            n = len(y[m][count1:count2])
            stderr[i].append(np.std(y[m][count1:count2])/np.sqrt(n))
            count1 = count2
                
    ind = np.arange(len(chosenchannels))
    width = 2
    _, ax = plt.subplots()
    fccolours = itertools.cycle(['red', 'gray','green', 'blue', 'cyan',
        'magenta','lawngreen','darkorange','maroon','orchid'])
        
    for i in range(len(paths)):
        paths[i] = os.path.basename(paths[i])[1:-5]
        test = list(paths[i])
        for z in range(len(test)):
            if test[z] == '_':
                test[z] = ' '
        paths[i] = ''.join(test)

    for i in range(len(paths)):
        fc = next(fccolours)
        ax.bar(ind + ind*(len(paths)*width) + width*i, means[i], width, 
            color=fc, alpha=.5, align='edge',
            error_kw=dict(ecolor='k', lw=1.5, capsize=3, capthick=1.5), 
            yerr=stderr[i], label=paths[i])#"mass scan {}".format(i)
        
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Mass/charge ratio")
    ax.set_xticks(ind + ind*(len(paths)*width) + width*(len(paths)/2))
    xlabels = [n[4:-2] for n in chosenchannels]
    ax.set_xticklabels(xlabels)
    ax.legend()
    plt.show()

def filenames_sort(filenames):
    entries = []
    for x in range(len(filenames)):
        entries.append(str(pd.read_excel(filenames[x], sheet_name="Time   Cycle")["Absolute Time"][0]) 
        + " " + filenames[x])
    entries.sort()
    filenames = [z.split()[2] for z in entries]
    return filenames

def plot_time_series():
    print("Plotting time series")

    #, constrained_layout=True)
    #ax1.grid(which='major', axis='both',color='skyblue',ls=':',lw=1)
    
    markercolours = itertools.cycle(['k','blue','r','magenta','midnightblue','darkorange'])
    linecolours = itertools.cycle(['grey','skyblue', 'orange','orchid', 'grey','g'])
    linestyles = itertools.cycle(['-'])
        
    all_channels, transmission = get_channels()
    channels = [i.get() for i in app.channels]
    filenames = filenames_sort(app.paths[0].get().split())
    ydata, ylabel, chosenchannels, rawsignal, filenames, ENvalues = get_ydata(channels, all_channels, transmission, filenames)
    
    for i in range(len(filenames)):
        filenames[i] = filenames_sort(filenames[i])

    if app.params[5].get() == "Normalised signal intensity":
        print("Plotting normalised signal intensity")

    xdata, absolute_time, xlabel, cycles_perxmins, date = get_xdata(app.params[6].get(), filenames)

    if app.yaxis2.get() == 1:
        if len(ydata) < 2:
            print("ydata needs to be length 2")
        else:
            print("correct ydata length")
            fig1, ax1 = plt.subplots(figsize=(20,10))
            ax2 = ax1.twinx()
            ax2.spines["right"].set_position(("axes", 1))
            ax1.set(xlabel=xlabel)
    
    else:
        fig1, ax1 = plt.subplots(figsize=(20,10))
        ax1.set(xlabel=xlabel, ylabel=ylabel)#, title=title)
    
    if app.EN.get() == 1:
#E/N == (Udrift*Vm*Td*p0)/(NA*T0*pd*ldrift)
#ENvalues = [Udrift, Tdrift, pdrift] everything in cgs!!
        Vm = 22414
        p0 = 1013.25
        NA = 6.022e23
        T0 = 273.15
        ldrift = 9.3
        print(ENvalues[0][0], ENvalues[1][0], ENvalues[2][0])
        EN = (np.array(ENvalues[0])*Vm*np.array(ENvalues[1])*p0)/(ldrift*NA*T0*np.array(ENvalues[2]))*1e17
        ysmooth = smooth(EN, cycles_perxmins[0])
##Need to make a second axis on the right hand side
        ax2 = ax1.twinx()
        ax2.set_ylabel("E/N (Td)")
        ax2.spines["right"].set_position(("axes", 1))
        ax1.plot(xdata[0][0], ydata[0][0], label= "E/N", lw=2, color="maroon")
        ax2.grid(False)
        ax2.plot(xdata[0], EN, '.',  ms=2, color='r', lw=1)
        ax2.plot(xdata[0][-1+cycles_perxmins[0]//2:-cycles_perxmins[0]//2], ysmooth,
        lw=2, color="maroon", label="E/N (Td)", linestyle='-')
    
    #PROBLEM: The cycles_per_xmins is per each file NOT each m/z - atm the cycles per x mins
    #is taken from the first file that's read in
    for index in range(len(chosenchannels)):
        #Could the ydata and xdata be called when the file is loaded in?
        #i.e. have the lists ready to go after load data is pressed?
        print("There are " + str(cycles_perxmins[index]) + " cycles per " 
              + str(app.params[1].get()) + " seconds.")
    
        mc = next(markercolours)
        lc = next(linecolours)
        ls = next(linestyles)

        compound = ','
        if chosenchannels[index] == "m/z 69":
            compound = r', isoprene (C$_5$H$_9^+$)'
        elif chosenchannels[index] == "m/z 87":
            compound = r', pentanal (C$_5$H$_{11}$O$^+$)'
        elif chosenchannels[index] == "m/z 59":
            compound = r', propanal (C$_3$H$_{7}$O$^+$)'
        elif chosenchannels[index] == "m/z 75":
            compound = r', diethyl-ether (C$_4$H$_{11}$O$^+$)'
        elif chosenchannels[index] == "m/z 43":
            compound = r', 1-propanol (C$_3$H$_{7}^+$)'
        elif chosenchannels[index] == "m/z 30":
            compound = r', NO$^+$'
        elif chosenchannels[index] == "m/z 95":
            compound = r', phenol (C$_6$H$_7$O$^+$)'
        elif chosenchannels[index] == "m/z 79":
            compound = r', benzene (C$_6$H$_7^+$)'

        series_label = (chosenchannels[index] + compound +'\n ' + str(app.params[1].get())
        + ' second moving average')

        ysmooth = smooth(ydata[index], cycles_perxmins[index])
        #Find the difference between the baseline correction and the minimum value of the baseline correction to apply to the 
        #data.
      #  ydata[index] = ydata[index] - min(ysmooth)
      #  ysmooth = ysmooth - min(ysmooth)
        
        if app.usebaseline.get() == 1:
            #Note: Fix this baseline so that the background signal is taken away. 
            #Need to know from the user when the blank is, then take this away from the ydata

            ##Use the times to get the correct index - similar to the readme.
            times = app.baselineparams[0].get().split('-')
            times = [n.strip() for n in times]
            
            datetime_object1 = datetime.datetime.strptime(
                date + ' ' + times[0], '%Y-%m-%d %H:%M:%S')
            cycle1 = bisect.bisect_right(absolute_time[index], datetime_object1)

            datetime_object2 = datetime.datetime.strptime(
                date + ' ' + times[1], '%Y-%m-%d %H:%M:%S')
            cycle2 = bisect.bisect_right(absolute_time[index], datetime_object2)
            baseline = ydata[index][cycle1:cycle2]
            baseline_mean = np.mean(baseline)
            print(baseline_mean)
            data = ydata[index] - baseline_mean
            factors = app.baselineparams[1].get().split(',')
            factors = [float(n) for n in factors]
            print(chosenchannels[index], factors[index])
            ydata[index] = np.array(data)/factors[index]

            ysmooth = smooth(ydata[index], cycles_perxmins[index])

            if app.logy.get() == 0:
                if app.yaxis2.get() == 1 and index==1:
                    ax1.plot(xdata[index][0], ydata[index][0], lw=2, color=lc, label=series_label,
                        linestyle=ls)
                    ax2.plot(xdata[index], ydata[index], '.',  ms=2, color=mc, lw=1)#,alpha=0.5,
                    ax2.plot(xdata[index][-1+cycles_perxmins[index]//2:-cycles_perxmins[index]//2], ysmooth, lw=2, color=lc, label=series_label,
                         linestyle=ls)
                    ax2.yaxis.set_minor_formatter(ScalarFormatter())
                    ax2.yaxis.set_major_formatter(ScalarFormatter())
                    ax2.set(ylabel=chosenchannels[index] + " " + ylabel )
                    ax1.set(ylabel=chosenchannels[index-1] + " " + ylabel)
                    ax1.grid(False)
                    ax2.grid(False)
                    if app.params[4].get() != "ymin - ymax":
                        ymin = float(app.params[4].get().split('-')[2])
                        ymax = float(app.params[4].get().split('-')[3])
                        ax2.set_ylim(ymin, ymax)
                else:
                    ax1.plot(xdata[index], ydata[index], '.',  ms=2, color=mc, lw=1)#,alpha=0.5,
                    ax1.plot(xdata[index][-1+cycles_perxmins[index]//2:-cycles_perxmins[index]//2], ysmooth, lw=2, color=lc, label=series_label,
                         linestyle=ls)
                    ax1.yaxis.set_minor_formatter(ScalarFormatter())
                    ax1.yaxis.set_major_formatter(ScalarFormatter())
                    if app.params[4].get() != "ymin - ymax":
                        ymin = float(app.params[4].get().split('-')[0])
                        ymax = float(app.params[4].get().split('-')[1])
                        ax1.set_ylim(ymin, ymax)

            elif app.logy.get() == 1:
                if app.yaxis2.get() == 1 and index==1:
                    ax1.plot(xdata[index][0], ydata[index][0], lw=2, color=lc, label=series_label,
                        linestyle=ls)
                    ax2.semilogy(xdata[index], ydata[index], '.',  ms=2, color=mc, lw=1)#,alpha=0.5,
                    ax2.semilogy(xdata[index][-1+cycles_perxmins[index]//2:-cycles_perxmins[index]//2], ysmooth, lw=2, color=lc, label=series_label,
                         linestyle=ls)
                    ax2.yaxis.set_minor_formatter(ScalarFormatter())
                    ax2.yaxis.set_major_formatter(ScalarFormatter())
                    ax2.set(ylabel=chosenchannels[index] + " " + ylabel )
                    ax1.set(ylabel=chosenchannels[index-1] + " " + ylabel)
                    ax1.grid(False)
                    ax2.grid(False)
                    if app.params[4].get() != "ymin - ymax":
                        ymin = float(app.params[4].get().split('-')[2])
                        ymax = float(app.params[4].get().split('-')[3])
                        ax2.set_ylim(ymin, ymax)
                else:
                    ax1.semilogy(xdata[index], ydata[index], '.',  ms=2, color=mc, lw=1)#,alpha=0.5,
                    ax1.semilogy(xdata[index][-1+cycles_perxmins[index]//2:-cycles_perxmins[index]//2], ysmooth, lw=2, color=lc, label=series_label,
                         linestyle=ls)
                    ax1.yaxis.set_minor_formatter(ScalarFormatter())
                    ax1.yaxis.set_major_formatter(ScalarFormatter())
                    if app.params[4].get() != "ymin - ymax":
                        ymin = float(app.params[4].get().split('-')[0])
                        ymax = float(app.params[4].get().split('-')[1])
                        ax1.set_ylim(ymin, ymax)
            
        elif app.usebaseline.get() == 0:
            if app.logy.get() == 0:
                if app.yaxis2.get() == 1 and index==1:
                    ax1.plot(xdata[index][0], ydata[index][0], lw=2, color=lc, label=series_label,
                        linestyle=ls)
                    ax2.plot(xdata[index], ydata[index], '.',  ms=2, color=mc, lw=1)#, alpha=0.5
                    ax2.plot(xdata[index][-1+cycles_perxmins[index]//2:-cycles_perxmins[index]//2], ysmooth, lw=2, color=lc, label=series_label,
                        linestyle=ls)
                    ax2.yaxis.set_minor_formatter(ScalarFormatter())
                    ax2.yaxis.set_major_formatter(ScalarFormatter())
                    ax2.set(ylabel=chosenchannels[index] + " " + ylabel )
                    ax1.set(ylabel=chosenchannels[index-1] + " " + ylabel)
                    ax1.grid(False)
                    ax2.grid(False)
                    if app.params[4].get() != "ymin - ymax":
                        ymin = float(app.params[4].get().split('-')[2])
                        ymax = float(app.params[4].get().split('-')[3])
                        ax2.set_ylim(ymin, ymax)
                else:
                    ax1.plot(xdata[index], ydata[index], '.',  ms=2, color=mc, lw=1)#, alpha=0.5
                    ax1.plot(xdata[index][-1+cycles_perxmins[index]//2:-cycles_perxmins[index]//2], ysmooth, lw=2, color=lc, label=series_label,
                        linestyle=ls)
                    ax1.yaxis.set_minor_formatter(ScalarFormatter())
                    ax1.yaxis.set_major_formatter(ScalarFormatter())
                    if app.params[4].get() != "ymin - ymax":
                        ymin = float(app.params[4].get().split('-')[0])
                        ymax = float(app.params[4].get().split('-')[1])
                        ax1.set_ylim(ymin, ymax)

            elif app.logy.get() == 1:
                if app.yaxis2.get() == 1 and index==1:
                    ax1.plot(xdata[index][0], ydata[index][0], lw=2, color=lc, label=series_label,
                        linestyle=ls)
                    ax2.semilogy(xdata[index], ydata[index], '.',  ms=2, color=mc, lw=1)#, alpha=0.5
                    ax2.semilogy(xdata[index][-1+cycles_perxmins[index]//2:-cycles_perxmins[index]//2], ysmooth, lw=2, color=lc, label=series_label,
                        linestyle=ls)
                    ax2.set(ylabel=chosenchannels[index] + " " + ylabel )
                    ax1.set(ylabel=chosenchannels[index-1] + " " + ylabel)
                    ax1.grid(False)
                    ax2.grid(False)
                    if app.params[4].get() != "ymin - ymax":
                        ymin = float(app.params[4].get().split('-')[2])
                        ymax = float(app.params[4].get().split('-')[3])
                        ax2.set_ylim(ymin, ymax)
                else:
                    ax1.semilogy(xdata[index], ydata[index], '.',  ms=2, color=mc, lw=1)#, alpha=0.5
                    ax1.semilogy(xdata[index][-1+cycles_perxmins[index]//2:-cycles_perxmins[index]//2], ysmooth, lw=2, color=lc, label=series_label,
                        linestyle=ls)
                    if app.params[4].get() != "ymin - ymax":
                        ymin = float(app.params[4].get().split('-')[0])
                        ymax = float(app.params[4].get().split('-')[1])
                        ax1.set_ylim(ymin, ymax)

    if app.paths[1].get() != '':
        use_readme(date, absolute_time[0], xdata[0], ydata, ax1, chosenchannels, rawsignal)
    
    if app.paths[2].get() != '' and app.paths[3].get() == '':
        #app.paths[2] == broadband
        plot_spectroscopy(app.paths[2].get(), app.params[2].get(), date, ax1, ["m", "b"], [1.2], 4)
        
    if app.paths[3].get() != '' and app.paths[2].get() == '':
        #app.paths[3] == OHN2
        plot_spectroscopy(app.paths[3].get(), app.params[3].get(), date, ax1, ["r", "g"], [1.2], 4)
    
    if app.paths[2].get() != '' and app.paths[3].get() != '':
        plot_spectroscopy(app.paths[3].get(), app.params[3].get(), date, ax1, ["r", "g"], [1.2], 4)
        plot_spectroscopy(app.paths[2].get(), app.params[2].get(), date, ax1, ["m", "b"], [1.4], 6)

   # leg = ax1.legend(ncol=3).get_frame() #
   # leg.set_alpha(0.8)
   # leg.set_edgecolor('white')
    
#    figManager = plt.get_current_fig_manager()
 #   figManager.app.showMaximized()
   # ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
   # plt.show()

   # fig1.savefig('myimage.svg', format='svg', dpi=1200)
   # fig1.savefig('myimage.eps', format='eps', dpi=1200)
   # fig2.savefig('myimage2.svg', format='svg', dpi=1200)
   # fig2.savefig('myimage2.eps', format='eps', dpi=1200)

    def ConctoDen(x):
#Accepts a concentration in ppb and returns the species density in cm^-3
        p = 101325
        k = 1.380649e-23
        T = 293.15
        n = x*1e-9*1e-6*p/(k*T)
        return n
    def DentoConc(n):
        #Accepts a species density in cm^-3 and returns the concentration in ppb
        p = 101325
        k = 1.380649e-23
        T = 293.15
        x = n*1e9*1e6*k*T/p
        return x

    title = date + '_experiment_' + chosenchannels[0].replace("/", "") 
    title.replace(" ","_")
    fig1.canvas.set_window_title(title)
    
    if app.specden.get() == 1:
        #ax1.set(ylabel=" m/z 95 Concentration (ppbv)")
        #ax2.set(ylabel="m/z 79 Concentration (ppbv)")
        secax = ax1.twinx()
        secax.spines["left"].set_position(("axes", -0.15))
        secax.yaxis.set_label_position('left')
        secax.yaxis.set_ticks_position('left')

#Need to calculate the most "useful" order of 10 for the offset calculation.
#Problem at the minute is using max - leads to nan issues.
#The nan issues are caused by the normalised cps calculations, the first few measurements
#of H3O+ (mz21) are zeroes, leading to divide by zero errors - meaning there is inf entries
# in the calculated normalised cps measurements.

        tmp = [x for x in ydata[0] if str(x) != 'nan']
        tmp = [x for x in tmp if float(x) != 0]
        tmp = np.asarray(tmp)
        tmp = tmp[tmp < 1e20]
            
        offset1 = "11" #str(int(np.log10(ConctoDen(tmp[0]))))
        
        offset = 1e11 #float("1e"+offset1)
        secax.set(ylabel=r"Species density ( $\times10^{" + offset1 + "}$cm$^{-3}$)")
        secax.plot(xdata[0], [ConctoDen(n)/offset for n in ydata[0]], color=None, alpha=0)

     #   scale_x = float( "1e" + offset1)
     #   print(scale_x)
     #   ticks_x = ticker.FuncFormatter(lambda x, pos: '{:.2f}'.format(x/scale_x))
     #   secax.yaxis.set_major_formatter(ticks_x)
        secax.grid(False)
        if app.params[4].get() != "ymin - ymax":
            ylimits = ax1.get_ylim()
#Need to get the values from the left hand side plot first.
            secax.set_ylim(ConctoDen(ylimits[0])/offset,ConctoDen(ylimits[1])/offset)

   # secax.grid(b=None)

    #secax.set(ylabel="Species density (cm$^{-3}$)")

   # secax = ax1.secondary_yaxis('right', functions=(ConctoDen,DentoConc))
   # ax1.figure.canvas.draw()
   # offset = secax.yaxis.get_major_formatter().get_offset()
   # secax.set(ylabel="Species density (" + offset + "cm$^{-3}$)")
   # secax.yaxis.offsetText.set_visible(False)

   # scale_x = float("1e"+offset[11:13])
   # ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:.3f}'.format(x/scale_x))
   # secax.yaxis.set_major_formatter(ticks_x)
    
    if app.params[6].get() == "Absolute Time":
        print("setting xaxis format to absolute time")
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))

  #  leg = ax1.legend(ncol=2).get_frame()
  #  leg.set_alpha(1)
  #  leg.set_edgecolor('None')
  #  fig1.tight_layout()
    ax1.legend(fancybox=False, ncol=2, framealpha=1, edgecolor="black", loc='lower left', bbox_to_anchor=(0, 1))#
    plt.show()

def plot_spectroscopy(path, lines, date, ax, colors, position, index):
#index = 4 for OH, N2 lines and 6 when Oxy atom lines
    filenames = os.listdir(path)
    filenames.sort()
    filenames = [path + '/' + f for f in filenames]
    lines = [l.strip() for l in lines.split(",")]
    xdata = []
    ydata = [[] for _ in range(len(lines))]
    
    for f in filenames:        
        with open(f) as file:
            data = file.readlines()
            for line in data:
                if "Date:" in line:
                    xdata.append(datetime.datetime.strptime(date + ' ' + line.strip().split()[4], "%Y-%m-%d %H:%M:%S"))

    for i in range(len(lines)):
        for f in filenames:
            with open(f) as file:
                data = file.readlines()
                for line in data:
                    if str(lines[i]) in line:
                        ydata[i].append(float(line.strip().split()[1]))    

    if app.params[6].get() == "Relative Time":
        reloffset = date + ' ' + app.params[0].get()
        reloffset = datetime.datetime.strptime(reloffset.strip(),
                                               '%Y-%m-%d %H:%M:%S')               
        xdata = [(item-reloffset).total_seconds()/60.0 for item in xdata]

    if len(lines) == 1:
        comp = ''
        if str(lines[0]) == "308.92":
            comp = ", OH"
        if str(lines[0]) == "844.66":
            comp = ", O"

        ax2 = ax.twinx()
        ax2.set_ylabel("Absolute Irradiance ($\mu$W/cm$^2$/nm)\n {} nm".format(str(lines[0])), color=colors[0])
        mean = np.mean(ydata[0][:3])
        ax2.plot(xdata, [x-mean for x in ydata[0]], label= "{}nm peak".format(str(lines[0])),lw=2,color=colors[0])
        ax2.spines["right"].set_position(("axes", position[0]))
        ax.plot(xdata[0], ydata[0][0], label= "{}~nm emission line".format(str(lines[0])) + comp,lw=2,color=colors[0])
        ax2.grid(False)
        if app.params[4].get() != "ymin - ymax":
            ymin = float(app.params[4].get().split('-')[index])
            ymax = float(app.params[4].get().split('-')[index+1])
            ax2.set_ylim(ymin, ymax)
        if app.params[6].get() == "Absolute Time":
            ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    
    elif len(lines) == 2:
        comp1 = ''
        comp2 = ''
        if str(lines[0]) == "308.92":
            comp1 = ", OH"
        if str(lines[1]) == "357.56":
            comp2 = r", N$_2$"
        if str(lines[1]) == "844.66":
            comp2 = ", O"

        ax2 = ax.twinx()
        ax2.set_ylabel("Absolute Irradiance ($\mu$W/cm$^2$/nm)\n {} nm".format(str(lines[0])), color=colors[0])
        mean = np.mean(ydata[0][:3])
        ax2.plot(xdata, [x-mean for x in ydata[0]], label= "{}nm peak".format(str(lines[0])),lw=2,color=colors[0])
        ax.plot(xdata[0], ydata[0][0], label= "{}~nm emission line".format(str(lines[0]))+comp1,lw=2,color=colors[0])
        ax2.grid(False)
        ax3 = ax.twinx()
        ax3.spines["right"].set_position(("axes", 1.2))
        ax3.set_ylabel("Absolute Irradiance ($\mu$W/cm$^2$/nm)\n {} nm".format(str(lines[1])), color=colors[1])
        mean = np.mean(ydata[1][:3])
        ax3.plot(xdata, [x-mean for x in ydata[1]], label= "{}nm peak".format(str(lines[1])),lw=2,color=colors[1])
        ax.plot(xdata[0], ydata[1][0], label= "{}~nm emission line".format(str(lines[1]))+comp2,lw=2,color=colors[1])
        ax3.grid(False)
        if app.params[6].get() == "Absolute Time":
            ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            ax3.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
        if app.params[4].get() != "ymin - ymax":
            ymin = float(app.params[4].get().split('-')[index])
            ymax = float(app.params[4].get().split('-')[index+1])
            ax2.set_ylim(ymin, ymax)
            ymin = float(app.params[4].get().split('-')[index+2])
            ymax = float(app.params[4].get().split('-')[index+3])
            ax3.set_ylim(ymin, ymax)
    
    elif len(lines) == 3:
        ax2 = ax.twinx()
        ax2.set_ylabel("Absolute Irradiance ($\mu$W/cm$^2$/nm)\n {} nm".format(str(lines[0])), color=colors[0])
        ax2.plot(xdata, ydata[0], label= "{}nm peak".format(str(lines[0])),lw=2,color=colors[0])
        ax.plot(xdata[0], ydata[1][0], label= "{}~nm emission line".format(str(lines[1])),lw=2,color=colors[0])
        ax2.grid(False)
        ax3 = ax.twinx()
        ax3.spines["right"].set_position(("axes", 1.2))
        ax3.set_ylabel("Absolute Irradiance ($\mu$W/cm$^2$/nm)\n {} nm".format(str(lines[1])), color=colors[1])
        ax3.plot(xdata, ydata[1], ls=':', label= "{}nm peak".format(str(lines[1])),lw=2,color=colors[1])
        ax.plot(xdata[0], ydata[1][0], label= "{}~nm emission line".format(str(lines[1])),lw=2,color=colors[1])
        ax3.grid(False)
        ax4 = ax.twinx()
        ax4.spines["right"].set_position(("axes", 1.3))
        ax4.set_ylabel("Absolute Irradiance ($\mu$W/cm$^2$/nm)\n {} nm".format(str(lines[2])), color=colors[2])
        ax4.plot(xdata, ydata[2], ls=':', label= "{}nm peak".format(str(lines[2])),lw=2,color=colors[2])
        ax.plot(xdata[0], ydata[1][0], label= "{}~nm emission line".format(str(lines[1])),lw=2,color=colors[2])
        ax4.grid(False)
        if app.params[6].get() == "Absolute Time":
            ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            ax3.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
            ax4.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))

def baseline_als(y, lam, p, niter=10):
  L = len(y)
  D = sparse.csc_matrix(np.diff(np.eye(L), 2))
  w = np.ones(L)
  for i in range(niter):
    W = sparse.spdiags(w, 0, L, L)
    Z = W + lam * D.dot(D.transpose())
    z = spsolve(Z, w*y)
    w = p * (y > z) + (1-p) * (y < z)
  return z

def trim_axs(axs, N):
    """little helper to massage the axs list to have correct length..."""
   # axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    return axs[:N]

def use_readme(date, absolute_time, xdata, ydata, ax, chosenchannels, rawsignal):

    indices = []
    times = []
    events = []
    cycles = []
    eventcycles = []
    y_calibdata = []
    y_errcalibdata = []
    rawsig_calibdata = []
    rawsig_stddev = []
    stddevs = []
    
    with open(app.paths[1].get()) as file:
        searchlines = file.readlines()
        for i, line in enumerate(searchlines):
            if "----" in line:
                indices.append(i)
                              
    with open(app.paths[1].get()) as file:
        searchlines = file.readlines()
        for l in searchlines[indices[2]+1:indices[3]]:
            time = l.split(',')
            time = [m.strip() for m in time]
            times.append(time)
      #  for l in searchlines[indices[4]+1:indices[5]]:
      #      ax.set(title = date + ' ' + l)      
      #  for l in searchlines[indices[6]+1:indices[7]]:
      #      dilution = l.split(',')
      #      dilution = [m.strip() for m in dilution]
      #  for l in searchlines[indices[4]+1:indices[5]]:
      #      event = l.split(',')
      #      event = [m.strip() for m in event]
      #      events.append(l)
            
   # ax.plot(absolute_time[0],[0],color='purple',label="OH spectra taken")
   # for x in range(len(events)):
   #     datetime_object1 = datetime.datetime.strptime(
   #         date + ' ' + events[x].strip(), '%Y-%m-%d %H:%M:%S')
   #     cycle = bisect.bisect_right(absolute_time, datetime_object1)
   #     eventcycles.append(cycle)
   # for x in range(len(eventcycles)):
   #     ax.axvline(xdata[eventcycles[x]],lw=1.5, color = 'purple')

    for x in range(len(times)):
        datetime_object1 = datetime.datetime.strptime(
            date + ' ' + times[x][0].strip(), '%Y-%m-%d %H:%M:%S')
        cycle1 = bisect.bisect_right(absolute_time, datetime_object1)

        datetime_object2 = datetime.datetime.strptime(
            date + ' ' + times[x][1].strip(), '%Y-%m-%d %H:%M:%S')
        cycle2 = bisect.bisect_right(absolute_time, datetime_object2)
#Part above is the normal action, part below is to restrict to only 100 measurement cycles.
       # cycle2 = cycle1 + 100
        cycles.append([cycle1, cycle2])
                                   
    cycle_labels = ['(' + string.ascii_uppercase[x] + ')' for x in range(len(cycles))]
      #  + '. ' + times[x][2] for x in range(len(cycles))]

       # for x in range(10):   
       #     time = datetime.datetime.strptime(reloffset, 
       #         '%Y-%m-%d %H:%M:%S') + datetime.timedelta(seconds=x*35)   
       #     cycle = bisect.bisect_right(absolute_time, time)
       #     ax1.axvline(xdata[cycle],lw=1, color = 'r')
    arrow_height=float(app.arrow_heightent.get())
    arrow_heightoffset=float(app.arrow_heightoffsetent.get())

    for x in range(len(cycles)):

        ax.axvline(xdata[(cycles[x])[0]],lw=1.5, color = times[x][3])#, label=cycle_labels[x])
        ax.axvline(xdata[(cycles[x])[1]],lw=1.5, color = times[x][3])
       # ax.plot(xdata[0],[0], label=cycle_labels[x], color='white')
            
        ax.annotate('', xy=(xdata[(cycles[x])[0]], arrow_height+x*arrow_heightoffset), 
            xytext=(xdata[(cycles[x])[1]], arrow_height+x*arrow_heightoffset), 
            arrowprops=dict(connectionstyle="arc3", arrowstyle="-", color=times[x][3]))
        
        ax.annotate(
            cycle_labels[x],((xdata[(cycles[x][0]+cycles[x][1])//2]), arrow_height + x*arrow_heightoffset + 0.5*arrow_heightoffset),
            horizontalalignment="center").draggable()
    
    dilution = []
    for x in range(len(cycles)):
        dilution.append(float(times[x][2]))

    for y in range(len(chosenchannels)):
        filename = (os.path.dirname(os.path.abspath(app.paths[0].get().split()[0])) + 
                   date + '_' + app.params[5].get() + '_' +
                   chosenchannels[y].replace("/", "") + ".txt")

        with open(filename, 'w') as fi:
            fi.write("#Label, xvalue, Mean, Standard deviation, Number of samples, Standard error of the mean, std err (%), Mean (cps)\n")
            for x in range(len(cycles)):
                tempydata =  ydata[y][cycles[x][0]:cycles[x][1]]
                temprawsig =  rawsignal[y][cycles[x][0]:cycles[x][1]]      
                mean = np.mean(tempydata)
                stddev = np.std(tempydata)
                rawsigmean = np.mean(temprawsig)
                rawsigdev = np.std(temprawsig)
                n = len(tempydata)
                stderr = stddev/np.sqrt(n)
                stderrper = stderr/mean*100
                print(string.ascii_uppercase[x] + '.\nmean: ' + str(mean) 
                    + '\nstd. dev.: ' + str(stddev) + '\nn: ' + str(n) 
                    + '\nstd. err.: ' + str(stderr) + '\nstd. err. %: ' + str(stderrper)
                    + '\nmean(cps): ' + str(rawsigmean) + '\n')
                fi.write(string.ascii_uppercase[x] + ', '+ str(dilution[x]) 
                    + ', ' + str(mean) 
                    + ', ' + str(stddev) + ', ' + str(n) 
                    + ', ' + str(stderr) + ', ' + str(stderrper) 
                    + ', ' + str(rawsigmean) + '\n')
###NOTE: part below for the change in samples & decrease in error
           #     for n in ydata[cycles[x][0]:cycles[x][1]]:
           #         fi.write("{0}\n".format(n))

                y_calibdata.append(mean)
                rawsig_calibdata.append(rawsigmean)
                rawsig_stddev.append(rawsigdev)
                y_errcalibdata.append(stderr)
                stddevs.append(stddev)

    if app.params[8].get() == 1:

        #Need to get a title, use polyfit for getting m, m_sigma, c and c_sigma
        #and use linregress for the r_value. Need also to get dilution/predicted
        #concentration from file        

        
        
       # filename = os.path.dirname(os.path.abspath(app.paths[0].get().split()[0])) + date + ' ' + chosenchannels[0].replace("/", "") + ".txt"
       # with open(filename, 'w') as fi:
       #     fi.write("# Label, dilution, Mean, Standard deviation, Number of points averaged over, Standard error of the mean (calculated using numpy)\n")
       #     for y in range(len(ydata)):
       #         for x in range(len(cycles)):         
       #             mean = np.mean(ydata[y][cycles[x][0]:cycles[x][1]])
       #             stddev = np.std(ydata[y][cycles[x][0]:cycles[x][1]])
       #             n = len(ydata[y][cycles[x][0]:cycles[x][1]])
       #             stderr = stddev/np.sqrt(n)                
       #             fi.write(string.ascii_uppercase[x] + ', ' 
       #                 + str(dilution[x]) + ', ' + str(mean) 
       #                 + ', ' + str(stddev) + ', ' + str(n) 
       #                 + ', ' + str(stderr) + ', ' + str(stderrper) + '\n')

        #dilution = 
        #[0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
        #[1.25e-3, 1.5e-3, 1.75e-3, 2e-3, 2.25e-3, 2.5e-3] # changed here
        #[0, 1.25e-3, 1.5e-3, 1.75e-3, 2e-3, 2.25e-3, 2.5e-3] for benzene
        #[0, 5e-3, 0.01, 0.015, 0.02, 0.025, 0.03] # correct for diethyl ether
        fig2, ax2 = plt.subplots(figsize=(15,15))#, constrained_layout=True)
        title = date + "_calibration" + "_" + chosenchannels[0].replace("/", "") 
        title.replace(" ","_")
        fig2.canvas.set_window_title(title)
        ax2.errorbar(dilution, y_calibdata, yerr=y_errcalibdata,
                        fmt='x',lw=1.5, ms=7, mew=1.5,capsize=5, 
                        color='k', capthick=1.5, label=chosenchannels[0])
        
        if len(dilution) == 4:
            rows, cols = 2, 2
        else:
            rows = 3
            if len(dilution)%rows == 0:
                cols = len(dilution)//rows
            else:
                cols = len(dilution)//rows + 1
        fig3 = plt.figure()
        ax = fig3.add_subplot(111) 
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
        ax.grid(False)
        co = 1
        axs = []
        for c in range(cols):
            for r in range(rows):
                c = c + 1
                r = r + 1
                axs.append(fig3.add_subplot(cols, rows, co))
                co += 1

     #   fig3, axs = plt.subplots(rows, cols, constrained_layout=True)
        title = date + "_probability_densities" + "_" + chosenchannels[0].replace("/", "") 
        title.replace(" ","_")
        fig3.canvas.set_window_title(title)
        
        if app.params[5].get() == "Normalised signal intensity":
            ax.set_xlabel("Normalised signal intensity (ncps)")
        elif app.params[5].get() == "Concentration":
            ax.set_xlabel("Measured concentration (ppbv)")

        ax.set_ylabel("Probability density")
        axs = trim_axs(axs, len(dilution))
        
        cos = range(len(dilution))
        for ax, co in zip(axs, cos):
            _, bins, _ = ax.hist(ydata[cycles[co][0]:cycles[co][1]], bins=30, density=True)#,label=cycle_labels[co])
            y = ((1 / (np.sqrt(2 * np.pi) * stddevs[co])) *  np.exp(-0.5 * (1 / stddevs[co] * (bins - y_calibdata[co]))**2))
            
            if co == 0:
                ax.axvline(y_calibdata[co], color='k', label=r"$\bar{x}$")
                ax.axvline(y_calibdata[co] + stddevs[co], color='r', label=r"$\sigma$")
                ax.plot(bins, y, '--', label=r"$f(x)$")
                ax.legend(fancybox=False, framealpha=1, edgecolor="black", loc='lower center', bbox_to_anchor=(0.5, 1), ncol=3)
                
            else:
                ax.plot(bins, y, '--')
                ax.axvline(y_calibdata[co], color='k')
                ax.axvline(y_calibdata[co] + stddevs[co], color='r')
            ax.axvline(y_calibdata[co] - stddevs[co], color='r')
            ax.annotate(s="{}".format(cycle_labels[co]), xy=(0.9, 0.8),xycoords="axes fraction")
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.0e'))
            ax.locator_params(nbins=5, axis='x')
            ax.tick_params(axis='y', which='both', labelleft=False)
        
        
        for n in range(len(dilution)):
            ax2.annotate(cycle_labels[n],(dilution[n],y_calibdata[n]+0.5)).draggable()

        slope, intercept, r_value, _, std_err = linregress(dilution, y_calibdata)
        f, V  = np.polyfit(dilution, y_calibdata, 1, cov=True, w=stddevs)
        title = date + ' characterisation'

#The caption for the calibration need to be "s=XX\pmYY \n i=XX\pmYY \n r^2=XX"
        ax2.plot(dilution, np.polyval(np.poly1d(f), dilution), 
            'r-', lw=1.5, label=(r"$s=$"+str(f[0])[:4]+"$\pm$"+str(np.sqrt(V[0][0]))[:3]+
            "\n"+"$i=23$"#+str(f[1])[:2]
            +"$\pm$"+str(np.sqrt(V[1][1]))[:1]+"\n"+"r$^2=$"+str(r_value)[:4]))
       
        ax2.axhspan(ymin=-50, ymax=y_calibdata[0]+3*stddevs[0], color="grey", alpha=0.5, label="LOD")
        
        #Linear regression does not take into account the standard deviation of each point.
       # ax2.plot(dilution, np.polyval([slope, intercept], dilution),
       #     'b--', lw=1.5, label='scipy linear regression')

        ax2.tick_params(which='both',direction='in',top=True, right=True)
      #  ax2.grid(which='major', axis='both',color='skyblue',ls=':',lw=1)
        if app.params[5].get() == "Normalised signal intensity":
            ax2.set_ylabel("Normalised signal intensity (ncps)")
        elif app.params[5].get() == "Concentration":
            ax2.set_ylabel("Measured concentration (ppbv)")
        #ax2.set(xlabel=r'Dilution factor ($\frac{\phi_{VOC}}{\phi{tot}}$)')
        ax2.set(xlabel="Concentration (ppbv)")
        ax2.legend(fancybox=False, framealpha=1, edgecolor="black")
        fig3.tight_layout()

        with open(filename, 'a') as fi:
            fi.write('# ' + chosenchannels[0] + '\n')
            fi.write('# y = mx + c\n')
            fi.write('# numpy polyfit\n')
            fi.write('# m = {}, sigma_m = {}\n'.format(f[0],np.sqrt(V[0][0])))
            fi.write('# c = {}, sigma_c = {}\n'.format(f[1],np.sqrt(V[1][1])))
            fi.write('# scipy linregress does not take standard deviation of each point into account.\n')
            fi.write('# scipy linregress\n')
            fi.write('# m = {}\n'.format(slope))
            fi.write('# c = {}\n'.format(intercept))
            fi.write('# r_value = {}\n'.format(r_value))
            fi.write('# sigma_m = {}\n'.format(std_err))
        
     #   fig4, ax4 = plt.subplots()
     #   print(y_calibdata, rawsig_calibdata)
     #   ax4.loglog(y_calibdata, rawsig_calibdata)

root = Tk()
app = Application(master=root)
app.mainloop()
root.destroy()