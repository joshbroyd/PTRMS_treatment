# Need to plot wavelength vs average y-value and plot error 
# bars/boxplot to see where the largest changes are.
import os
import matplotlib.pyplot as plt
import numpy as np


folder_name = "/home/jgb509/Documents/CRM/Spectroscopy/20190404_OHN2_DE_cl2/plasma_on/"
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

plt.show()

                