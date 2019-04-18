from tkinter import (Frame, Tk, Button, filedialog, Entry, OptionMenu, 
                    StringVar, Label, IntVar, Checkbutton)
import pandas as pd

class Application(Frame):
    
    def load_data(self):
        self.tickboxes = []
        self.channel_entries = [[] for _ in range(3)] 
        channels = get_channels()
        for m in range(1, 3):
            for n in range(len(channels[m])):            
                entry = IntVar()
                tickbox = Checkbutton(self, text=channels[m][n],
                            variable=entry)
                tickbox.grid(row=6+n, column=m-1, sticky="NESW") #m+(len(channel_keys[0])/14 + 1)
                self.channel_entries[m].append(entry)
                self.tickboxes.append(tickbox)
        
        if len(channels[0]) > 19:
            mz_channels = [n.split() for n in channels[0]]
            available_mz = ("m/z " + str(mz_channels[0][1]) + " to " 
                + str(mz_channels[-1][1]) + " are available")
            usr_inst = ("Input individual comma separated values\n"
                "or a range, hypen separated.")
            Label(self, text=available_mz).grid(row=5, column=2, sticky="NESW")
            Label(self, text=usr_inst).grid(row=7, column=2, sticky="NESW")
            self.channel_entries[0] = Entry(self)   
            self.channel_entries[0].insert("0", str(mz_channels[0][1]) + "-" 
                + str(mz_channels[-1][1]))
            self.channel_entries[0].grid(row=6, column=2)
        
        elif len(channels[0]) < 20:
            for n in range(len(channels[0])):
                x = n//14 + 2
                y = n + 6 - (n//14 * 14)           
                entry = IntVar()
                tickbox = Checkbutton(self, text=channels[0][n], 
                        variable=entry)
                tickbox.grid(row=y, column=x, sticky="NESW")
                self.channel_entries[0].append(entry)
                self.tickboxes.append(tickbox)

    def createWidgets(self):
        Button(self, text="Excel files", command=self.get_excel_paths).grid(
            row=0, column=0, sticky="NSEW")
        self.file_entry = Entry(self)
        self.file_entry.grid(row=0, column=1, sticky="NSEW")
        
        Label(self, text="y-axis format").grid(row=0, column=2, sticky="NSEW")
        self.yformat = StringVar(self)
        self.yformat.set("Concentration")
        OptionMenu(self, self.yformat, *["Concentration", 
            "Raw signal intensities"]).grid(row=0, column=3, sticky="NSEW")

        Label(self, text="x-axis format").grid(row=1, column=2, sticky="NSEW")
        self.xformat = StringVar(self)
        self.xformat.set("Relative Time")
        OptionMenu(self, self.xformat, *["Cycle number", "Absolute Time", 
            "Relative Time"]).grid(row=1, column=3, sticky="NSEW")

        Button(self, text="Readme file", command=get_readme_path).grid(
            row=1, column=0, sticky="NSEW")
        self.readme_path = Entry(self)
        self.readme_path.grid(row=1, column=1, sticky="NSEW")

        Button(self, text="Spectroscopy folder", 
            command=get_spec_path).grid(row=2, column=0, sticky="NSEW")
        self.spec_path = Entry(self)
        self.spec_path.grid(row=2, column=1, sticky="NSEW")

        Button(self, text="Load data", 
            command=self.load_data).grid(row=3, column=0, sticky="NSEW")

        Button(self, text="Plot time series", 
            command=plot_time_series).grid(row=4, column=0, sticky="NSEW")

        Label(self, text="Moving average\nduration (seconds)").grid(row=4,
            column=2, sticky="NESW")
        self.mov_av = Entry(self)
        self.mov_av.insert("0", 120.0)
        self.mov_av.grid(row=4, column=3, sticky="NSEW")

        Label(self, text="Absolute time (hh:mm) \nwhen relative time = 0").grid(
            row=3, column=2, sticky="NESW")
        self.rel_offset = Entry(self)
        self.rel_offset.grid(row=3, column=3, sticky="NSEW")

        self.calib = IntVar()
        Checkbutton(self, text="Calibration", 
            variable=self.calib).grid(column=2,row=2, sticky="NESW")

        Button(self, text="Plot mass scan", 
            command=plot_mass_scan).grid(row=4, column=1, sticky="NSEW")

        Button(self, text="QUIT", fg="red", command=self.quit).grid(
            row=5, column=0, sticky="NSEW")
            
    def get_excel_paths(self):
        #This needs to be here to insert the pathnames to the entry, 
        #does it though? if the entry itself is an attribute to this class?
        self.file_paths = filedialog.askopenfilenames(
            initialdir="/home/jgb509/Documents/CRM/Data")
        self.file_entry.insert("0", self.file_paths)
        print(self.file_paths)

    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.createWidgets()



def get_readme_path():
    app.readme_path.insert("0", filedialog.askopenfilename(
        initialdir="/home/jgb509/Documents/CRM/Data"))

def get_spec_path():
    app.spec_path.insert("0", filedialog.askdirectory(
        initialdir="/home/jgb509/Documents/CRM/Data"))

def get_channels():
#Accepts a list of filenames and returns the full list of channel names (keys) 
#from the first excel file that's read in

    all_channels = []
    sheetnames = ["Raw signal intensities", "Instrument",
                  "Reaction conditions"]
    paths = app.file_paths
    print(paths)    
    for sheetname in sheetnames:
        data = pd.read_excel(paths[0], sheet_name=sheetname)
        channels = [str(x) for x in list(data.keys())]
        all_channels.append(channels)
    
    #Need to check the channels of the other files in turn as well.

    return all_channels

def plot_time_series():
    print("plotting time series")

def plot_mass_scan():
    print("plotting mass scan")

root = Tk()
app = Application(master=root)
app.mainloop()
root.destroy()