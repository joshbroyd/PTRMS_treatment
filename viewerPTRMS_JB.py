from tkinter import (Tk, Frame, Button, filedialog, Entry, IntVar, Checkbutton,
Label, StringVar, OptionMenu)
import matplotlib as mpl
mpl.use('qt5agg')
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import sys, itertools, os, datetime, bisect, string
from scipy.stats import linregress
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#Changes the font and fontsize of the graphs 
if __name__ == "__main__":

    fontsize = 20
    params = {'backend':'qt5agg',
            'text.latex.preamble':['\\usepackage{gensymb}'],
            'axes.labelsize':fontsize,
            'axes.titlesize':fontsize,
            'font.size':fontsize,
            'legend.fontsize':fontsize-5,
            'xtick.labelsize':fontsize-5,
            'ytick.labelsize':fontsize-5,
            'font.family':'serif'}
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
        # "[4] = Graph title
        # "[5] = y-axis format
        # "[6] = x-axis format
        # "[7] = time series or mass scan
        # "[8] = calibration
        r, c = [3, 4, 1, 3, 4], [3, 3, 5, 5, 5]
        labeltext = ["Absolute time (hh:mm:ss) \nwhen trel=0",
                     "Moving average\nduration (seconds)",
                     "Broadband lines \n (O peak at 777.25nm and 844.66nm):",
                     ("OHN2 lines (OH peak at 308.92nm, \nN2 peaks at " 
                      "336.30nm, 357.56nm):"), "Graph title:"]
        defaults = ['', 120.0, "777.25, 844.66", "308.92, 336.30, 357.56", 
                    "experiment"]
        for i in range(len(labeltext)):
            Label(self, text=labeltext[i]).grid(column=c[i]-1, row=r[i])
            self.params[i] = Entry(self)
            self.params[i].grid(column=c[i],row=r[i])
            self.params[i].insert("0", defaults[i])
        
        menutext = ["y-axis format", "x-axis format",""]
        defaultoption = ["Concentration", "Absolute Time", "Time series"]
        options = [["Concentration", "Raw signal intensities"], 
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
        Button(self, text="Load data", command=self.loaddata).grid(
            row=2, column=0, sticky="NSEW")
        Button(self, text="Plot", command=plot).grid(
            row=3,column=0, sticky="NSEW")
        Button(self, text="Exit", command=self.quit).grid(
            row=4, column=0, sticky="NSEW")

    def getdatapath(self):
        path = list(filedialog.askopenfilenames(
            title="Choose data files",
            initialdir="/home/jgb509/Documents/CRM/Data"))
        self.paths[0].insert("0", path)

    def getreadmepath(self):
        path = filedialog.askopenfilename(
            title="Choose readme file",
            initialdir="/home/jgb509/Documents/CRM/Data")
        self.paths[1].insert("0", path)

    def getbroadpath(self):
        path = filedialog.askdirectory(
            title='Choose broadband spectroscopy folder', 
            initialdir="/home/jgb509/Documents/CRM/Spectroscopy")
        self.paths[2].insert("0", path)    

    def getOHN2path(self):
        path = filedialog.askdirectory(
            title='Choose OHN2 spectroscopy folder', 
            initialdir="/home/jgb509/Documents/CRM/Spectroscopy")
        self.paths[3].insert("0", path)

    def loaddata(self):
        all_channels = get_channels()
        self.tickboxes = []
        self.channels = [[] for _ in range(len(all_channels))]
        #self.channels[0] = user chosen m/z channels
        # "[1] = Recorded instrument settings
        # "[2] = Measured instrument conditions
        for n in range(1, 3):
            for m in range(len(all_channels[n])):
                entry = IntVar()
                tickbox = Checkbutton(self, text=all_channels[n][m], 
                                      variable=entry)
                tickbox.grid(row=5+m, column=n-1, sticky="NESW")
                self.channels[n].append(entry)
                self.tickboxes.append(tickbox)
        
        if self.params[7].get() == "Time series":
            for n in range(len(all_channels[0])):
                x = n//14 + 2
                y = n + 6 - (n//14 * 14)           
                entry = IntVar()
                tickbox = Checkbutton(self, text=all_channels[0][n], 
                                      variable=entry)
                tickbox.grid(row=y, column=x, sticky="NESW")
                self.channels[0].append(entry)
                self.tickboxes.append(tickbox)
        
        elif self.params[7].get() == "Mass scan":
            mz_channels = [n.split() for n in all_channels[0]]
            available_mz = ("m/z " + str(mz_channels[0][1]) + " to " 
                + str(mz_channels[-1][1]) + " are available")
            usr_inst = ("Input individual comma separated values\n"
                        "or a range, hypen separated.")
            Label(self, text=available_mz).grid(row=5, column=2, sticky="NESW")
            Label(self, text=usr_inst).grid(row=7, column=2, sticky="NESW")
            self.channels[0] = Entry(self)   
            self.channels[0].insert("0", str(mz_channels[0][1]) + "-" 
                + str(mz_channels[-1][1]))
            self.channels[0].grid(row=6, column=2)

def get_channels():
#Accepts a list of filenames and returns the full list of channel names (keys) 
#from the first excel file that's read in

    all_channels = []
    sheetnames = ["Raw signal intensities", "Instrument",
                  "Reaction conditions"]
    path = app.paths[0].get().split()          
    for sheetname in sheetnames:
        data = pd.read_excel(path[0], sheet_name=sheetname)
        channels = [str(x) for x in list(data.keys())]
        all_channels.append(channels)
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

def get_xdata(xoption, reloffset):
#Returns xdata list and xlabel depending on xoption.
    
    if reloffset != " " and xoption == "Relative Time":
        reloffset = datetime.datetime.strptime(reloffset.strip(),
                                               '%Y-%m-%d %H:%M:%S')
    
    if reloffset == " " and xoption == "Relative Time":
        sys.exit("Need to choose a starting time for relative time")
        
    x = []
    for f in app.paths[0].get().split():
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

def get_ydata(yoption, channels, channel_keys):
#Args are corresponding to raw counts or concentration and the list of channel 
#keys. Returns a list of lists of ydata from the corresponding channels and 
#yoption and the yaxis label. 

    chosenchannels = []
    for n in range(len(channels)):
        if channels[n] == 1:
            chosenchannels.append(channel_keys[n])
 
    y = [[] for _ in range(len(chosenchannels))]      
    for f in app.paths[0].get().split():
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
        
    return y, ylabel, chosenchannels

def smooth(y, box_pts):
#Smooth function to return the moving average of length box_pts from list y.     
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def plot():
    if app.params[7].get() == "Time series":
        plot_time_series()
        for tickbox in app.tickboxes:
            tickbox.destroy()

    elif app.params[7].get() == "Mass scan":
        plot_mass_scan()
        for tickbox in app.tickboxes:
            tickbox.destroy()

def plot_mass_scan():
    print("plotting mass scan")
                
    all_channels = get_channels()
    channels = np.zeros(len(all_channels[0]), dtype=int)
    n_channels = convertparam(app.channels[0].get())      
    n_channels = ["m/z " + str(i) for i in n_channels]
    for n in n_channels:
        channels[all_channels[0].index(n)] = 1
               
    y, ylabel, chosenchannels = get_ydata(app.params[5].get(), channels, 
        all_channels[0])
    print(ylabel, chosenchannels)
        
    lengths = []
    paths = app.paths[0].get()
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
        paths[i] = os.path.basename(paths[i])[:-5]
        
    for i in range(len(paths)):
        fc = next(fccolours)
        ax.bar(ind + ind*(len(paths)*width) + width*i, means[i], width, 
            color=fc, alpha=.5, align='edge',
            error_kw=dict(ecolor='k', lw=1.5, capsize=3, capthick=1.5), 
            yerr=stderr[i], label=paths[i])
        
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Mass/charge ratio")
    ax.set_xticks(ind + ind*(len(paths)*width) + width*(len(paths)/2))
    xlabels = [n[4:-2] for n in chosenchannels]
    ax.set_xticklabels(xlabels)
    ax.legend()

    plt.show()

def plot_time_series():
    print("plotting time series")

    _, ax1 = plt.subplots(figsize=(20,10))
    if app.params[6].get() == "Absolute Time":
        print("setting xaxis format to absolute time")
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))                  
                
    markercolours = itertools.cycle(['k','lightgreen','r','magenta',
                                     'midnightblue','darkorange'])
    linecolours = itertools.cycle(['grey','g','maroon','orchid', 'skyblue',
                                   'orange'])
    linestyles = itertools.cycle(['-', '--', ':', '-.'])
        
    all_channels = get_channels()
    channels = [i.get() for i in app.channels[0]]

    ydata, ylabel, chosenchannels = get_ydata(app.params[5].get(), channels, 
        all_channels[0])
    absolute_time = get_xdata("Absolute Time", app.params[0].get())[0]
    date = datetime.datetime.strftime(absolute_time[0], '%Y-%m-%d')
    reloffset = date + ' ' + app.params[0].get()
    xdata, xlabel = get_xdata(app.params[6].get(), reloffset)   
    rel_time = get_xdata("Relative Time", 
        datetime.datetime.strftime(absolute_time[0], '%Y-%m-%d %H:%M:%S'))[0]
    cycles_perxmins = int(np.argmin([abs(element*60 - float(app.params[1].get()))
                          for element in rel_time]))                          
    print("There are " + str(cycles_perxmins) + " cycles per " 
          + str(app.params[1].get()) + " seconds.")

    for index in range(len(ydata)):
        mc = next(markercolours)
        lc = next(linecolours)
        ls = next(linestyles)
        series_label = (chosenchannels[index] + ',\n ' + str(app.params[1].get())
        + ' second moving average')
        ysmooth = smooth(ydata[index], cycles_perxmins)
        ax1.plot(xdata, ydata[index], 'o',  ms=2, color=mc)
        ax1.plot(xdata, ysmooth, lw=2, color=lc, label=series_label,
                 linestyle=ls)

    title = date + ' ' + app.params[4].get()
    ax1.set(xlabel=xlabel, ylabel=ylabel, title=title)
    if app.paths[1].get() != '':
        use_readme(date, absolute_time, xdata, ydata, ax1, chosenchannels)
    
    if app.paths[2].get() != '':
        plot_spectroscopy(app.paths[2].get(), app.params[2].get(), date, ax1)
    
    if app.paths[3].get() != '':
        plot_spectroscopy(app.paths[3].get(), app.params[3].get(), date, ax1)
    
    ax1.grid(which='major', axis='both',color='skyblue',ls=':',lw=1)
    ax1.yaxis.set_minor_formatter(ScalarFormatter())
    ax1.yaxis.set_major_formatter(ScalarFormatter())
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

    leg = ax1.legend(ncol=3).get_frame()
    leg.set_alpha(0.8)
    leg.set_edgecolor('white')
    plt.show()

def plot_spectroscopy(path, lines, date, ax):

    filenames = os.listdir(path)
    filenames.sort()
    filenames = [path + '/' + f for f in filenames] 
    xdata = []
    ydata = [[] for _ in range(len(lines))]
    lines = convertparam(lines)

    for f in filenames:        
        with open(f) as file:
            data = file.readlines()
            for line in data:
                if "Date:" in line:
                    xdata.append(datetime.datetime.strptime(date + ' ' + line.strip().split()[4], "%Y-%m-%d %H:%M:%S"))
                for i in range(len(lines)):
                    if str(lines[i]) in line:
                        ydata[i].append(float(line.strip().split()[1]))    

    if app.params[6].get() == "Relative Time":
        reloffset = date + ' ' + app.params[0].get()
        reloffset = datetime.datetime.strptime(reloffset.strip(),
                                               '%Y-%m-%d %H:%M:%S')               
        xdata = [(item-reloffset).total_seconds()/60.0 for item in xdata]

  #  print(app.params[6].get(), reloffset)

    if len(lines) == 1:
        ax2 = ax.twinx()
        ax2.set_ylabel("Absolute Irradiance ($\mu$W/cm$^2$/nm)\n {} nm".format(str(lines[0])), color='r')
        ax2.plot(xdata, ydata[0], label= "{}nm peak".format(str(lines[0])),lw=2,color='r')
        ax.plot(xdata[0], ydata[0][0], label= "{}nm emission line".format(str(lines[0])),lw=2,color='r')
    
    elif len(lines) == 2:
        ax2 = ax.twinx()
        ax2.set_ylabel("Absolute Irradiance ($\mu$W/cm$^2$/nm)\n {} nm".format(str(lines[0])), color='r')
        ax2.plot(xdata, ydata[0], label= "{}nm peak".format(str(lines[0])),lw=2,color='r')
        ax3 = ax.twinx()
        ax3.spines["right"].set_position(("axes", 1.1))
        ax3.set_ylabel("Absolute Irradiance ($\mu$W/cm$^2$/nm)\n {} nm".format(str(lines[1])), color='g')
        ax3.plot(xdata, ydata[1], ls=':', label= "{}nm peak".format(str(lines[1])),lw=2,color='g')


def use_readme(date, absolute_time, xdata, ydata, ax, chosenchannels):

    indices = []
    times = []
    events = []
    cycles = []
    eventcycles = []
    y_calibdata = []
    y_errcalibdata = []
    stddevs = []
    fig2, ax2 = plt.subplots(figsize=(15,15))
    
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
        for l in searchlines[indices[4]+1:indices[5]]:
            ax.set(title = date + ' ' + l)      
        for l in searchlines[indices[6]+1:indices[7]]:
            dilution = l.split(',')
            dilution = [m.strip() for m in dilution]
        for l in searchlines[indices[8]+1:indices[9]]:
          #  event = l.split(',')
          #  event = [m.strip() for m in event]
            events.append(l)
            
  #  ax1.plot((params[11])[0],[0],color='purple',label="N$_2$ temperature taken")
  #  for x in range(len(events)):
  #      datetime_object1 = datetime.datetime.strptime(
  #          params[9] + ' ' + events[x].strip(), '%Y-%m-%d %H:%M:%S')
  #      cycle = bisect.bisect_right(params[10], datetime_object1)
  #      eventcycles.append(cycle)
  #   for x in range(len(eventcycles)):
  #      ax1.axvline(params[11][eventcycles[x]],lw=1.5, color = 'purple')

    for x in range(len(times)):
        datetime_object1 = datetime.datetime.strptime(
            date + ' ' + times[x][0].strip(), '%Y-%m-%d %H:%M:%S')
        cycle1 = bisect.bisect_right(absolute_time, datetime_object1)

        datetime_object2 = datetime.datetime.strptime(
            date + ' ' + times[x][1].strip(), '%Y-%m-%d %H:%M:%S')
        cycle2 = bisect.bisect_right(absolute_time, datetime_object2)
        cycles.append([cycle1, cycle2])
                                   
    cycle_labels = [string.ascii_uppercase[x] 
        + '. ' + times[x][2] for x in range(len(cycles))]

       # for x in range(10):   
       #     time = datetime.datetime.strptime(reloffset, 
       #         '%Y-%m-%d %H:%M:%S') + datetime.timedelta(seconds=x*35)   
       #     cycle = bisect.bisect_right(absolute_time, time)
       #     ax1.axvline(xdata[cycle],lw=1, color = 'r')

    
    
    for x in range(len(cycles)):
        #can change this  len(cycles) to the labels you want to show
                 
        arrow_height=sum(ydata[0])/len(ydata[0]) - x + 2
             
        ax.axvline(xdata[(cycles[x])[0]],lw=1.5, color = times[x][3], label=cycle_labels[x])
        ax.axvline(xdata[(cycles[x])[1]],lw=1.5, color = times[x][3])
    #    ax1.plot([0],[0], label=cycle_labels[x], color='white')
            
        ax.annotate('', xy=(xdata[(cycles[x])[0]], arrow_height-0.1), 
            xytext=(xdata[(cycles[x])[1]], arrow_height-0.1), 
            arrowprops=dict(connectionstyle="arc3", arrowstyle="<->", color=times[x][3]))
            
        ax.annotate(
            cycle_labels[x][0],((xdata[cycles[x][0]]), arrow_height))#+ 20
                                     
    for y in range(len(ydata)):
        for x in range(len(cycles)):         
            mean = np.mean(ydata[y][cycles[x][0]:cycles[x][1]])
            stddev = np.std(ydata[y][cycles[x][0]:cycles[x][1]])
            n = len(ydata[y][cycles[x][0]:cycles[x][1]])
            stderr = stddev/np.sqrt(n)                
            print(string.ascii_uppercase[x] + '.\nmean: ' + str(mean) 
                + '\nstd. dev.: ' + str(stddev) + '\nn: ' + str(n) 
                + '\nstd. err.: ' + str(stderr) + '\n')
            y_calibdata.append(mean)
            y_errcalibdata.append(stderr)
            stddevs.append(stddev)
    
    if app.params[8].get() == 1:

        #Need to get a title, use polyfit for getting m, m_sigma, c and c_sigma
        #and use linregress for the r_value. Need also to get dilution/predicted
        #concentration from file        

        dilution = []
        for x in range(len(cycles)):
            dilution.append(float(times[x][2]))
            
        #dilution = 
        #[0, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
        #[1.25e-3, 1.5e-3, 1.75e-3, 2e-3, 2.25e-3, 2.5e-3] # changed here
        #[0, 1.25e-3, 1.5e-3, 1.75e-3, 2e-3, 2.25e-3, 2.5e-3] for benzene
        #[0, 5e-3, 0.01, 0.015, 0.02, 0.025, 0.03] # correct for diethyl ether

        ax2.errorbar(dilution, y_calibdata, yerr=y_errcalibdata,
                        fmt='x',lw=1.5, ms=7, mew=1.5,capsize=5, 
                        color='k', capthick=1.5, label=chosenchannels[0])

        f, V  = np.polyfit(dilution, y_calibdata, 1, cov=True, w=stddevs)
        
        title = date + ' characterisation'
        ax2.plot(dilution, np.polyval(np.poly1d(f), dilution), 
            'r-', lw=1.5, label='Linear fit')

        ax2.tick_params(which='both',direction='in',top=True, right=True)
        ax2.grid(which='major', axis='both',color='skyblue',ls=':',lw=1)
        ax2.set(xlabel='Dilution factor', 
            ylabel = 'Measured\nconcentration (ppb)', title=title)        
        leg = ax2.legend().get_frame()
        leg.set_alpha(1)
        leg.set_edgecolor('white')

        slope, intercept, r_value, _, std_err = linregress(dilution, y_calibdata)

        filename = date + ' ' + chosenchannels[0].replace("/", "") + ".txt"
        with open(filename, 'w') as fi:
            fi.write(chosenchannels[0])
            fi.write('\ny = mx + c\n')
            fi.write('numpy polyfit\n')
            fi.write('m = {}, sigma_m = {}\n'.format(f[0],np.sqrt(V[0][0])))
            fi.write('c = {}, sigma_c = {}\n'.format(f[1],np.sqrt(V[1][1])))
            fi.write('scipy linregress\n')
            fi.write('m = {}\n'.format(slope))
            fi.write('c = {}\n'.format(intercept))
            fi.write('r_value = {}\n'.format(r_value))
            fi.write('sigma_m = {}\n'.format(std_err))

root = Tk()
app = Application(master=root)
app.mainloop()
root.destroy()
