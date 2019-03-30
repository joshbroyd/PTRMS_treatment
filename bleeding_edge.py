from tkinter import (Tk, Entry, filedialog, Button, Label, StringVar,
                     OptionMenu, IntVar,Checkbutton)
import os
import pandas as pd

def select_params():

    window_1 = Tk()

    #Entry for the data directory path
    e1 = Entry(window_1)
    e1.grid(row=0, column=1, sticky="ew")

    #OptionMenu for showing whether to plot MID or SCAN
    e2 = StringVar(window_1)
    e2.set("Time series")
    OptionMenu(window_1, e2, *["Time series", 
        "Mass scan"]).grid(column=2,row=0, sticky="W")
    
    #OptionMenu for showing whether to plot concentration or raw signal
    e3 = StringVar(window_1)
    e3.set("Concentration")
    OptionMenu(window_1, e3, *["Concentration", 
        "Raw signal intensities"]).grid(column=4,row=0, sticky="W")

 #The below need to be included only after choosing time series
    #OptionMenu for showing whether to plot against cycle number, absolute 
    #time or relative time.
  #  Label(window_1, text="x-axis format").grid(
  #      row=1, sticky="W")
  #  ent2 = StringVar(window_1)
  #  ent2.set("Absolute Time")
  #  OptionMenu(window_1, ent2, *["Cycle number", "Absolute Time",
  #      "Relative Time"]).grid(column=1,row=1, sticky="W")
  #Entry for the number of minutes to use for the moving average
  #  Label(window_1, text="Moving average\nduration (seconds)").grid(row=2,
  #      sticky="W")
  #  ent3 = Entry(window_1)
  #  ent3.insert("0", 120.0)
  #  ent3.grid(row=2, column=1, sticky="ew")

    def getdatadir():
        #Should be able to make a list of files.
        global datadir
        filez = filedialog.askopenfilenames(
            initialdir = "/home/jgb509/Documents/CRM/Data")
        print(filez)

      #  datadir = filedialog.askdirectory(
      #      initialdir = "/home/jgb509/Documents/CRM/Data")
       # e1.insert("0", datadir)


    def close():
        window_1.quit()
        window_1.destroy()

    Button(window_1, text="Data directory", command=getdatadir).grid(row=0, 
        sticky="ew")
    Button(window_1, text="Exit", command=close).grid(row=4, sticky="ew")


    window_1.mainloop()

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

    # Add in check to throw error if not all files have the same/ same order
    # of channels.

    #Should throw all available channels from the selected files, so that if
    # one file has one extra or one less channel its not a problem, the 
    #program should be able to find where the channel is.

    return all_channels

select_params()