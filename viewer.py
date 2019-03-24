import os, string, itertools, bisect, peakutils, sys
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('WXAgg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
from matplotlib.pyplot import cm
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from tkinter import (filedialog, Tk, IntVar, Checkbutton, Label, Entry, 
                     StringVar, OptionMenu, Button, mainloop, Text)

#Fix so that it no longer crashes when "exit" is pressed
#Option for dilution and titles?
#Bring back the automated peak detection?
#Make the GUI as close as possible to the actual PTRMS viewer?

#Changes the font and fontsize of the graphs 
if __name__ == "__main__":

    fontsize = 35
    params = {'backend':'WXAgg',
            'text.latex.preamble':['\\usepackage{gensymb}'],
            'axes.labelsize':fontsize,
            'axes.titlesize':fontsize,
            'font.size':fontsize,
            'legend.fontsize':fontsize-5,
            'xtick.labelsize':fontsize-5,
            'ytick.labelsize':fontsize-5,
            'font.family':'serif'}
    mpl.rcParams.update(params)

def select_param():

    window_1 = Tk()
    
    def getchannelsts():     
        global all_channels, channel_keys, tickboxes
        tickboxes = []
        channel_keys = get_channels()      
        all_channels = [[] for _ in range(3)]
        
        for n in range(len(channel_keys[0])):
            x = n//14 + 3
            y = n + 5 - (n//14 * 14)           
            entry = IntVar()
            tickbox = Checkbutton(window_1, text=channel_keys[0][n], 
                        variable=entry)
            tickbox.grid(row=y, column=x, sticky="W")
            all_channels[0].append(entry)
            tickboxes.append(tickbox)
           
        for m in range(1, 3):
            for n in range(len(channel_keys[m])):            
                entry = IntVar()
                tickbox = Checkbutton(window_1, text=channel_keys[m][n],
                            variable=entry)
                tickbox.grid(row=5+n, column=m, 
                            sticky="W") #m+(len(channel_keys[0])/14 + 1)
                all_channels[m].append(entry)
                tickboxes.append(tickbox)
        
    def getdatadir():
        global datadir
        datadir = filedialog.askdirectory(
            initialdir = "/home/jgb509/Documents/CRM/PTR-MS/Data")
        ent5.insert("0", datadir)
        
    def getreadmedir():
        global readmedir
        readmedir = filedialog.askopenfilename(
            initialdir = "/home/jgb509/Documents/CRM/PTR-MS/Data")
        ent4.insert("0", readmedir)
    
    def timeseriesgetvals():       
        global all_channels
        reloffset = ent8.get()   
        yoption = ent1.get()
        yoption = [yoption, "Instrument", "Reaction conditions"]        
        xoption = ent2.get()
        secs_to_av = float(ent3.get())
        readmedir = ent4.get()
        caloption = int(ent6.get())
        
        for n in range(3):
            channels = [int(ent.get()) for ent in all_channels[n]]
            if any(channels) == 1:
                plotting(channel_keys[n], channels, yoption[0], xoption,
                     secs_to_av, readmedir, yoption[n], caloption, reloffset)

        for tickbox in tickboxes:
            tickbox.destroy()
    
    def getchannelsms():
        global ent7
        all_channels = get_channels()
        all_channels = [n.split() for n in all_channels[0]]
        available_mz = ("m/z " + str(all_channels[0][1]) + " to " 
            + str(all_channels[-1][1]) + " are available")
        usr_inst = ("Input individual comma separated values\n"
            "or a range, hypen separated.")
        Label(window_1, text=available_mz).grid(row=5, column=1)
        Label(window_1, text=usr_inst).grid(row=7, column=1)
        ent7 = Entry(window_1)   
        ent7.insert("0", str(all_channels[0][1]) + "-" 
            + str(all_channels[-1][1]))
        ent7.grid(row=6, column=1)
    
    def massscangetvals():
        global ent7
                
        yoption = ent1.get()
        all_channels = get_channels()
        channels = np.zeros(len(all_channels[0]), dtype=int)
        n_channels = convertparam(ent7.get())      
        n_channels = ["m/z " + str(i) for i in n_channels]
        for n in n_channels:
            channels[all_channels[0].index(n)] = 1
               
        y, ylabel, chosenchannels = get_ydata(yoption, channels, 
            all_channels[0], yoption)
   
        print(ylabel, chosenchannels)
        
        files = file_retrieval()
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
                
        print(means)        
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
        
        
    def close():
        window_1.quit()
        window_1.destroy()

    #OptionMenu for showing whether to plot concentration or raw signal
    Label(window_1, text="y-axis format").grid(row=0, sticky="W")
    ent1 = StringVar(window_1)
    ent1.set("Concentration")
    OptionMenu(window_1, ent1, *["Concentration", 
        "Raw signal intensities"]).grid(column=1,row=0, sticky="W")
    
    #OptionMenu for showing whether to plot against cycle number, absolute 
    #time or relative time.
    Label(window_1, text="x-axis format").grid(
        row=1, sticky="W")
    ent2 = StringVar(window_1)
    ent2.set("Absolute Time")
    OptionMenu(window_1, ent2, *["Cycle number", "Absolute Time",
        "Relative Time"]).grid(column=1,row=1, sticky="W")
    
    #Entry for the number of minutes to use for the moving average
    Label(window_1, text="Moving average\nduration (seconds)").grid(row=2,
        sticky="W")
    ent3 = Entry(window_1)
    ent3.insert("0", 120.0)
    ent3.grid(row=2, column=1, sticky="ew")
        
    #Entry for the readme file path
    ent4 = Entry(window_1)
    ent4.grid(row=4,column=1,sticky="ew")
    
    #Entry for the data directory path
    ent5 = Entry(window_1)
    ent5.grid(row=3,column=1,sticky="ew")
        
    #Checkbox for calibration plot
    ent6 = IntVar()
    Checkbutton(window_1, text="Calibration", 
        variable=ent6).grid(column=2,row=2, sticky="W")
    
    #Entry for the relative time t = 0 point
    Label(window_1, text="Absolute time that the\nrelative time starts"
        " from (hh:mm:ss)").grid(row=0, column=2, sticky="W")
    ent8 = Entry(window_1)
    ent8.insert("0", "11:49:22")
    ent8.grid(row=1, column=2, sticky="W")
    
    #Buttons
    Button(window_1, text="Data directory", command=getdatadir).grid(row=3, 
        sticky="ew")
    Button(window_1, text="README path", command=getreadmedir).grid(row=4, 
        sticky="ew")
    Button(window_1, text="Load time series \n channel list", 
           command=getchannelsts).grid(row=5, sticky="ew")
    Button(window_1, text="Load mass scan \n channel list",
           command=getchannelsms).grid(row=6, sticky="ew")
    Button(window_1, text="Plot time series",
           command=timeseriesgetvals).grid(row=8, sticky="ew")
    Button(window_1, text="Plot mass scan",
           command=massscangetvals).grid(row=9, sticky="ew")
    Button(window_1, text="Exit", command=close).grid(row=10, sticky="ew")
    
    mainloop()

def file_retrieval():
#Returns list of filenames in order.

    files = os.listdir(datadir)
    files.sort()     
    files = [datadir + "/" + f for f in files]      
    return files
  
def get_channels():
#Accepts a list of filenames and returns the full list of channel names (keys) 
#from the first excel file that's read in

    files = file_retrieval()
    all_channels = []
    sheetnames = ["Raw signal intensities", "Instrument",
                  "Reaction conditions"]                 
    for sheetname in sheetnames:
        data = pd.read_excel(files[0], sheet_name=sheetname)
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

def get_ydata(yoption, channels, channel_keys, sheetname):
#Args are corresponding to raw counts or concentration and the list of channel 
#keys. Returns a list of lists of ydata from the corresponding channels and 
#yoption and the yaxis label. 

    files = file_retrieval()
    
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

    files = file_retrieval()
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

def plotting(channel_keys, channels, yoption, xoption,
                 secs_to_av, readmedir, sheetname, caloption, reloffset):
    
 #   _, ax1 = plt.subplots()
    
    if xoption == "Absolute Time":
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))                    
                
    markercolours = itertools.cycle(['k','lightgreen','r','magenta',
                                     'midnightblue','darkorange'])
    linecolours = itertools.cycle(['grey','g','maroon','orchid', 'skyblue',
                                   'orange'])
    linestyles = itertools.cycle(['-', '--', ':', '-.'])
    
    print(sheetname)
    
    ydata, ylabel, chosenchannels = get_ydata(yoption, channels, channel_keys, 
                                              sheetname)
    absolute_time = get_xdata("Absolute Time", reloffset)[0]
    date = datetime.datetime.strftime(absolute_time[0], '%Y-%m-%d')
    
    reloffset = date + ' ' + reloffset
    
    xdata, xlabel = get_xdata(xoption, reloffset)   
    rel_time = get_xdata("Relative Time", 
        datetime.datetime.strftime(absolute_time[0], '%Y-%m-%d %H:%M:%S'))[0]
    
    
    cycles_perxmins = int(np.argmin([abs(element*60 - secs_to_av)
                          for element in rel_time]))                          
    print("There are " + str(cycles_perxmins) + " cycles per " 
          + str(secs_to_av) + " seconds.")

    for index in range(len(ydata)):
        mc = next(markercolours)
        lc = next(linecolours)
        ls = next(linestyles)
        series_label = (chosenchannels[index] + ', ' + str(secs_to_av)
        + ' second moving average')
        ysmooth = smooth(ydata[index], cycles_perxmins)
        ax1.plot(xdata, ydata[index], 'o',  ms=2, color=mc)
        ax1.plot(xdata, ysmooth, lw=2, color=lc, label=series_label,
                 linestyle=ls)  
                 
  #  title = date + ' TMO canister measurement'
  #  plt.title(title)
  
    '''
    ROIT = ["2019-02-27 13:15:00", "2019-02-27 13:24:00", "2019-02-27 14:15:00", "2019-02-27 14:30:00", "2019-02-27 14:51:00", "2019-02-27 15:14:00","2019-02-27 13:22:00"]
    ROIL = ["Benzene on","Plasma on","Plasma off", "Benzene off", "Plasma on", "Plasma off"]
    ROIC = ["lawngreen", 'orange', 'orange', "lawngreen", 'orange', 'orange', 'b']
    ROI_cycles = []
    
    for x in range(len(ROIT)):   
        time = datetime.datetime.strptime(ROIT[x], '%Y-%m-%d %H:%M:%S')
        cycle = bisect.bisect_right(absolute_time, time)
        ROI_cycles.append(cycle)
    #    ax1.axvline(xdata[cycle],lw=2, color=ROIC[x] )
    
    print(ROI_cycles)
    
    arrow_height = 1500 #was 10 for the phenol channel
    
    ax1.annotate('', xy=(xdata[ROI_cycles[0]], arrow_height-500), #was arrow_height-5 for phenol channel
                xytext=(xdata[ROI_cycles[3]], arrow_height-500), 
                arrowprops=dict(arrowstyle="<->", color='limegreen',lw=2)) #benzene
    
    ax1.annotate('', xy=(xdata[ROI_cycles[1]], arrow_height), 
                xytext=(xdata[ROI_cycles[2]], arrow_height), 
                arrowprops=dict(arrowstyle="<->", color='darkorange',lw=2)) # plasma
    
    ax1.annotate('', xy=(xdata[ROI_cycles[4]], arrow_height), 
                xytext=(xdata[ROI_cycles[5]], arrow_height), 
                arrowprops=dict(arrowstyle="<->", color='darkorange',lw=2)) # plasma
                
    ax1.annotate('', xy=(xdata[ROI_cycles[6]], arrow_height+500), 
                xytext=(xdata[ROI_cycles[3]], arrow_height+500), 
                arrowprops=dict(arrowstyle="<->", color='b',lw=2)) # argon for benzene channel
    
    ax1.plot(xdata[0], [0], label= "Plasma 40W$_f$",lw=2,color='darkorange')
    ax1.plot(xdata[0], [0], label= "Benzene on",lw=2,color='limegreen')
    ax1.plot(xdata[0], [0], label= "Argon on",lw=2,color='b')'''        
    
    if readmedir != '':

        indices = []
        times = []
        cycles = []
        y_calibdata = []
        y_errcalibdata = []
        inverse_stddev = []

        README = readmedir
    
        with open(README) as file:
            searchlines = file.readlines()
            for i, line in enumerate(searchlines):
                if "----" in line:
                    indices.append(i)
                              
        with open(README) as file:
             searchlines = file.readlines()
             for l in searchlines[indices[2]+1:indices[3]]:
                 time = l.split(',')
                 time = [m.strip() for m in time]
                 times.append(time)

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
                 
            arrow_height=sum(ydata[0])/len(ydata[0])
             
            ax1.axvline(xdata[(cycles[x])[0]],lw=1.5, color = 'b')
            ax1.axvline(xdata[(cycles[x])[1]],lw=1.5, color = 'b')
          #  ax1.plot(xdata[0],[0], label=cycle_labels[x], color='white')
            
            ax1.annotate('', xy=(xdata[(cycles[x])[0]], arrow_height-0.1), 
                xytext=(xdata[(cycles[x])[1]], arrow_height-0.1), 
                arrowprops=dict(arrowstyle="<->", color='b'))
            
            ax1.annotate(
                cycle_labels[x][0],((xdata[cycles[x][0]
                    + 20]), 20))
            
                            
        for y in range(len(ydata)):
            print(chosenchannels[y])
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
                inverse_stddev.append(1/stddev)
    
        if caloption == 1:
            
            _, ax2 = plt.subplots()
            dilution = []
           # dilution = [1.25e-3, 1.5e-3, 1.75e-3, 2e-3, 2.25e-3, 2.5e-3] # changed here
            
            #[0, 1.25e-3, 1.5e-3, 1.75e-3, 2e-3, 2.25e-3, 2.5e-3] for benzene
            #[0, 5e-3, 0.01, 0.015, 0.02, 0.025, 0.03] # correct for diethyl ether
            for x in range(len(cycles)):
                dilution.append(float(times[x][2]))
            
            print(len(dilution), len(y_calibdata))         
            ax2.errorbar(dilution, y_calibdata, yerr=y_errcalibdata,
                                fmt='x',lw=1.5, ms=7, mew=1.5,capsize=5, 
                                color='k', label=chosenchannels[0], capthick=1.5)
        
            fitted_data = np.polyfit(dilution, y_calibdata, 1, full=True, 
                w=inverse_stddev)
        
            print(fitted_data)
        
            label = str(np.poly1d(fitted_data[0]))
             
            ax2.plot(dilution, np.polyval(np.poly1d(fitted_data[0]), dilution), 
                'r-', lw=1.5, label='y = ' + label[2:])

            ax2.tick_params(which='both',direction='in',top=True, right=True)
            ax2.set( xlabel='Dilution factor',
                    ylabel = 'Measured concentration (ppb)')#title=title,
                    
            leg = ax2.legend()
            leg.get_frame().set_alpha(1)
            leg.get_frame().set_edgecolor('white')
            ax2.grid(which='major', axis='both',color='skyblue',ls=':',lw=1)  
                                     
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.grid(which='major', axis='both',color='skyblue',ls=':',lw=1)
    ax1.yaxis.set_minor_formatter(ScalarFormatter())
    ax1.yaxis.set_major_formatter(ScalarFormatter())
    leg = ax1.legend() #ncol=2
    leg.get_frame().set_alpha(1)
    leg.get_frame().set_edgecolor('white')    
    plt.show()



#####
#Main code
#####

_, ax1 = plt.subplots()
'''folder_name = "1slm ar 40Wf with benzene/"

filenames = os.listdir(folder_name)
filenames.sort()
filenames = [folder_name + f for f in filenames]

# OH peak at 308.92nm, NH peak at 336.30nm, N2 peak at 357.56nm
 
xdata = []
ydata = [[] for _ in range(4)]

for f in filenames:        

    with open(f) as file:
        data = file.readlines()
        for line in data:
            if "Date:" in line:
                xdata.append(datetime.datetime.strptime("2019-02-27 " + line.strip().split()[4], "%Y-%m-%d %H:%M:%S"))
            if "308.92" in line:
                ydata[0].append(float(line.strip().split()[1]))
            if "336.30" in line:
                ydata[1].append(float(line.strip().split()[1]))
            if "357.56" in line:
                ydata[2].append(float(line.strip().split()[1]))
                

ax2 = ax1.twinx()
ax3 = ax1.twinx()
ax4 = ax1.twinx()

ax2.axes.get_yaxis().set_ticks([])
ax3.axes.get_yaxis().set_ticks([])
ax4.axes.get_yaxis().set_ticks([])

ax1.plot(xdata[0], [0], label= "308.92nm peak",lw=2,color='r')
ax1.plot(xdata[0], [0], ls='--', label= "336.30nm peak",lw=2,color='b')
ax1.plot(xdata[0], [0], ls=':', label= "357.56nm peak",lw=2,color='g')

ax2.plot(xdata, ydata[0], label= "308.92nm peak",lw=2,color='r')
ax3.plot(xdata, ydata[1], ls='--', label= "336.30nm peak",lw=2,color='b') #[:100] for the increase after just plasma on
ax4.plot(xdata, ydata[2], ls=':', label= "357.56nm peak",lw=2,color='g')


folder_name = "1slm ar 40Wf without benzene/"

filenames = os.listdir(folder_name)
filenames.sort()
filenames = [folder_name + f for f in filenames]

# OH peak at 308.92nm, NH peak at 336.30nm, N2 peak at 357.56nm
 
xdata = []
ydata = [[] for _ in range(4)]

for f in filenames:        

    with open(f) as file:
        data = file.readlines()
        for line in data:
            if "Date:" in line:
                xdata.append(datetime.datetime.strptime("2019-02-27 " + line.strip().split()[4], "%Y-%m-%d %H:%M:%S"))
            if "308.92" in line:
                ydata[0].append(float(line.strip().split()[1]))
            if "336.30" in line:
                ydata[1].append(float(line.strip().split()[1]))
            if "357.56" in line:
                ydata[2].append(float(line.strip().split()[1]))

ax2.plot(xdata, ydata[0],lw=2,color='r')
ax3.plot(xdata, ydata[1], ls='--',lw=2,color='b') #[:100] for the increase after just plasma on
ax4.plot(xdata, ydata[2], ls=':',lw=2,color='g')

#ax1.xlabel("Wavelength (nm)")
#ax1.ylabel("Absolute Irradiance ($\mu$W/cm$^{2}$/nm)")'''

select_param()




