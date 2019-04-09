
from tkinter import (filedialog, Tk, Button, Label, StringVar, OptionMenu, 
    Entry, IntVar, Checkbutton)
import matplotlib as mpl
mpl.use('qt5agg') #WXAgg
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import numpy as np
import sys, itertools, os, datetime, bisect, string
import pandas as pd
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from scipy.stats import linregress

#Make params into a dict.

#Changes the font and fontsize of the graphs 
if __name__ == "__main__":

    fontsize = 35
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

def select_param():
    
    global params, specdir
    params = [[] for _ in range(9)]

    print(len(params))

    commands = [getdatapaths, getreadmepath, loaddata, plot, close]
    commandtext = ["Add Excel Files", "Add Project File" , "Load data", "Plot",
        "Exit"] 
    for i in range(len(commandtext)):
        Button(window, text=commandtext[i], command=commands[i]).grid(row=i, 
        sticky="NESW")

    Button(window, text="Add Spectra Folder", command=getspecfolder).grid(row=5, 
        sticky="NESW")
    specdir = Entry(window)
    specdir.grid(column=1, row=5)


    Label(window, text="Moving average\nduration (seconds)").grid(row=4,
        column=2, sticky="NESW")
    Label(window, text="Absolute time (hh:mm) \nwhen relative time = 0").grid(
        row=3, column=2, sticky="NESW")

    paramsindex = range(3,7)
    rowindex = [3, 1, 0, 4]
    columnindex = [3, 1, 1, 3]
    for i in range(len(paramsindex)):
        params[paramsindex[i]] = Entry(window)
        params[paramsindex[i]].grid(column=columnindex[i], row=rowindex[i])

    menutext = ["y-axis format", "x-axis format",""]
    defaultoption = ["Concentration","Absolute Time","Time series"]
    options = [["Concentration", "Raw signal intensities"], ["Cycle number", 
                "Absolute Time", "Relative Time"],["Time series", "Mass scan"]]
    for i in range(3):
        Label(window, text=menutext[i]).grid(row=i, column=2, sticky="NESW")
        params[i] = StringVar(window)
        params[i].set(defaultoption[i])
        OptionMenu(window, params[i], *options[i]).grid(row=i, column=3, sticky="NESW")
            
    params[6].insert("0", 120.0)
    params[7] = IntVar()
    Checkbutton(window, text="Calibration", 
        variable=params[7]).grid(column=2,row=2, sticky="NESW")

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

def getspecfolder():
    path = filedialog.askdirectory(parent=window,
        title='Choose spectroscopy folder', 
        initialdir="/home/jgb509/Documents/CRM/Spectroscopy")
    specdir.insert("0", path)

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
                tickbox.grid(row=6+n, column=m-1, sticky="NESW") #m+(len(channel_keys[0])/14 + 1)
                params[8][m].append(entry)
                tickboxes.append(tickbox)
    
    if params[2].get() == "Time series":
        for n in range(len(channels[0])):
            x = n//14 + 2
            y = n + 6 - (n//14 * 14)           
            entry = IntVar()
            tickbox = Checkbutton(window, text=channels[0][n], 
                        variable=entry)
            tickbox.grid(row=y, column=x, sticky="NESW")
            params[8][0].append(entry)
            tickboxes.append(tickbox)

    elif params[2].get() == "Mass scan":
        mz_channels = [n.split() for n in channels[0]]
        available_mz = ("m/z " + str(mz_channels[0][1]) + " to " 
            + str(mz_channels[-1][1]) + " are available")
        usr_inst = ("Input individual comma separated values\n"
            "or a range, hypen separated.")
        Label(window, text=available_mz).grid(row=5, column=2, sticky="NESW")
        Label(window, text=usr_inst).grid(row=7, column=2, sticky="NESW")
        params[8][0] = Entry(window)   
        params[8][0].insert("0", str(mz_channels[0][1]) + "-" 
            + str(mz_channels[-1][1]))
        params[8][0].grid(row=6, column=2)

def get_channels():
#Accepts a list of filenames and returns the full list of channel names (keys) 
#from the first excel file that's read in

    all_channels = []
    sheetnames = ["Raw signal intensities", "Instrument",
                  "Reaction conditions"]                 
    for sheetname in sheetnames:
        data = pd.read_excel(paths[0], sheet_name=sheetname)
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

def plot():
    if params[2].get() == "Time series":
        plot_time_series()
        if specdir.get() != '':
            plot_spectroscopy()
        plt.show()

        for tickbox in tickboxes:
            tickbox.destroy()

    elif params[2].get() == "Mass scan":
        plot_mass_scan()

def plot_spectroscopy():

    folder_name = specdir.get() + '/'
   # print(folder_name)

    filenames = os.listdir(folder_name)
    filenames.sort()
    filenames = [folder_name + f for f in filenames]
  #  print(filenames)

# OH peak at 308.92nm, NH peak at 336.30nm, N2 peak at 357.56nm
 
    xdata = []
    ydata = [[] for _ in range(4)]

    for f in filenames:        

        with open(f) as file:
            data = file.readlines()
            for line in data:
                if "Date:" in line:
                    xdata.append(datetime.datetime.strptime(params[9] + ' ' + line.strip().split()[4], "%Y-%m-%d %H:%M:%S"))
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

'''
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



def plot_time_series():
    print("plotting time series")
    global ax1, ax2, chosenchannels

    fig1, ax1 = plt.subplots(figsize=(20,10))
    fig2, ax2 = plt.subplots(figsize=(15,15))

    if params[1].get() == "Absolute Time":
        print("setting xaxis format to absolute time")
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))                    
                
    markercolours = itertools.cycle(['k','lightgreen','r','magenta',
                                     'midnightblue','darkorange'])
    linecolours = itertools.cycle(['grey','g','maroon','orchid', 'skyblue',
                                   'orange'])
    linestyles = itertools.cycle(['-', '--', ':', '-.'])
        
    all_channels = get_channels()
    channels = [i.get() for i in params[8][0]]

    ydata, ylabel, chosenchannels = get_ydata(params[0].get(), channels, 
        all_channels[0])
    absolute_time = get_xdata("Absolute Time", params[3].get())[0]
    date = datetime.datetime.strftime(absolute_time[0], '%Y-%m-%d')
    
    reloffset = date + ' ' + params[3].get()
    
    xdata, xlabel = get_xdata(params[1].get(), reloffset)   
    rel_time = get_xdata("Relative Time", 
        datetime.datetime.strftime(absolute_time[0], '%Y-%m-%d %H:%M:%S'))[0]
     
    cycles_perxmins = int(np.argmin([abs(element*60 - float(params[6].get()))
                          for element in rel_time]))                          
    print("There are " + str(cycles_perxmins) + " cycles per " 
          + str(params[6].get()) + " seconds.")

    for index in range(len(ydata)):
        mc = next(markercolours)
        lc = next(linecolours)
        ls = next(linestyles)
        series_label = (chosenchannels[index] + ',\n ' + str(params[6].get())
        + ' second moving average')
        ysmooth = smooth(ydata[index], cycles_perxmins)
        ax1.plot(xdata, ydata[index], 'o',  ms=2, color=mc)
        ax1.plot(xdata, ysmooth, lw=2, color=lc, label=series_label,
                 linestyle=ls)

    title = date + ' experiment'
    ax1.set(xlabel=xlabel, ylabel=ylabel, title = title)
    if params[4].get() != '' or specdir != '':
        # DO NOT use append to append to params - the shape of params should be 
        # left alone
        print('adding things to params')
        params.append(date)
        params.append(absolute_time)
        params.append(xdata)
        params.append(ydata)
        if params[4].get() != '':
            use_readme()
    
    ax1.grid(which='major', axis='both',color='skyblue',ls=':',lw=1)
    ax1.yaxis.set_minor_formatter(ScalarFormatter())
    ax1.yaxis.set_major_formatter(ScalarFormatter())
    leg = ax1.legend().get_frame() #ncol=2
    leg.set_alpha(0.8)
    leg.set_edgecolor('white')
    
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
   # plt.show()

   # fig1.savefig('myimage.svg', format='svg', dpi=1200)
   # fig1.savefig('myimage.eps', format='eps', dpi=1200)
   # fig2.savefig('myimage2.svg', format='svg', dpi=1200)
   # fig2.savefig('myimage2.eps', format='eps', dpi=1200)
    
    

def plot_mass_scan():
    print("plotting mass scan")
                
    all_channels = get_channels()
    channels = np.zeros(len(all_channels[0]), dtype=int)
    n_channels = convertparam(params[8][0].get())      
    n_channels = ["m/z " + str(i) for i in n_channels]
    for n in n_channels:
        channels[all_channels[0].index(n)] = 1
               
    y, ylabel, chosenchannels = get_ydata(params[0].get(), channels, 
        all_channels[0])
    print(ylabel, chosenchannels)
        
    lengths = []
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

def use_readme():

    indices = []
    times = []
    cycles = []
    y_calibdata = []
    y_errcalibdata = []
    stddevs = []
    
    with open(params[4].get()) as file:
        searchlines = file.readlines()
        for i, line in enumerate(searchlines):
            if "----" in line:
                indices.append(i)
                              
    with open(params[4].get()) as file:
        searchlines = file.readlines()
        for l in searchlines[indices[2]+1:indices[3]]:
            time = l.split(',')
            time = [m.strip() for m in time]
            times.append(time)
        for l in searchlines[indices[4]+1:indices[5]]:
            ax1.set(title = params[9] + ' ' + l)

    for x in range(len(times)):
        datetime_object1 = datetime.datetime.strptime(
            params[9] + ' ' + times[x][0].strip(), '%Y-%m-%d %H:%M:%S')
        cycle1 = bisect.bisect_right(params[10], datetime_object1)

        datetime_object2 = datetime.datetime.strptime(
            params[9] + ' ' + times[x][1].strip(), '%Y-%m-%d %H:%M:%S')
        cycle2 = bisect.bisect_right(params[10], datetime_object2)
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
                 
        arrow_height=sum(params[12][0])/len(params[12][0])
             
        ax1.axvline(params[11][(cycles[x])[0]],lw=1.5, color = 'b')
        ax1.axvline(params[11][(cycles[x])[1]],lw=1.5, color = 'b')
          #  ax1.plot(xdata[0],[0], label=cycle_labels[x], color='white')
            
        ax1.annotate('', xy=(params[11][(cycles[x])[0]], arrow_height-0.1), 
            xytext=(params[11][(cycles[x])[1]], arrow_height-0.1), 
            arrowprops=dict(connectionstyle="arc3", arrowstyle="-", color='b'))
            
        ax1.annotate(
            cycle_labels[x][0],((params[11][cycles[x][0]
                + 20]), 20))
                                     
    for y in range(len(params[12])):
        for x in range(len(cycles)):         
            mean = np.mean(params[12][y][cycles[x][0]:cycles[x][1]])
            stddev = np.std(params[12][y][cycles[x][0]:cycles[x][1]])
            n = len(params[12][y][cycles[x][0]:cycles[x][1]])
            stderr = stddev/np.sqrt(n)                
            print(string.ascii_uppercase[x] + '.\nmean: ' + str(mean) 
                + '\nstd. dev.: ' + str(stddev) + '\nn: ' + str(n) 
                + '\nstd. err.: ' + str(stderr) + '\n')
            y_calibdata.append(mean)
            y_errcalibdata.append(stderr)
            stddevs.append(stddev)
    
    if params[7].get() == 1:

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
        
        title = params[9] + ' characterisation'
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

        filename = params[9] + ' ' + chosenchannels[0].replace("/", "") + ".txt"
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

        

window = Tk()
select_param()
window.mainloop()