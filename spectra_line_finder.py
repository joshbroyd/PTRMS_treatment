# Need to plot wavelength vs average y-value and plot error 
# bars/boxplot to see where the largest changes are.
import os
import matplotlib as mpl
mpl.use('qt5agg') #WXAgg
import matplotlib.pyplot as plt
import numpy as np



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


plt.rc('text', usetex=True)


#folder_name = "/home/jgb509/Documents/CRM/Spectroscopy/20190404_OHN2_DE_cl2/plasma_on/"
#folder_name = "/home/jgb509/Documents/CRM/Spectroscopy/20181207_broadband_readings/500sccm_argon_41Wf_1Wr/"
folder_name = "/home/jgb509/Documents/CRM/Spectroscopy/20190513_broadband_and_OHN2_O2_test/Broadband/Plasma_on/"
   # print(folder_name)

filenames = os.listdir(folder_name)
filenames.sort()
filenames = [folder_name + f for f in filenames]
  #  print(filenames)

#need a list of wavelengths and a list of lists for 
#each wavelength from each file
 
xdata = []
ydata = []

for f in [filenames[0]]:  
    with open(f) as file:
        data = file.readlines()
        for l in data[17:3665]:
            z = l.split()
            xdata.append(float(z[0]))

ydata = [[] for _ in range(len(xdata))]

for f in filenames:
    with open(f) as file:
        data = file.readlines()
        i = 0        
        for l in data[17:3665]:
            z = l.split()
            ydata[i].append(float(z[1]))
            i += 1

ydatastddev = [np.std(i) for i in ydata]
ydatameans = [np.mean(i) for i in ydata]
            
plt.errorbar(xdata, ydatameans, yerr=ydatastddev,
                        fmt='x',lw=1.5, ms=7, mew=1.5,capsize=5, 
                        color='k', capthick=1.5)
plt.xlabel("Wavelength (nm)")
plt.ylabel(r"Absolute Irradiance ($\mu$ W/cm$^2$/nm)") 
plt.show()

                
