import tkinter as tk
from tkinter import filedialog, ttk

#Add in ability to plot on more than one axis for time series, i.e. be
#able to plot m/z 21 and 75 on the same figure but overlapped, with
#their corresponding y-axis.

class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.grid()
        self.create_widgets()

    def create_widgets(self):

        tabControl = ttk.Notebook(root)          # Create Tab Control
        tab1 = ttk.Frame(tabControl)
        tab2 = ttk.Frame(tabControl)            # Create a tab 
        tabControl.add(tab1, text='Time Series')
        tabControl.add(tab2, text='Mass Scan')     # Add the tab
        tabControl.grid(row=3, column=0, columnspan=2)  # Pack to make visible

        menubar = tk.Menu(self.master)
        self.master.config(menu=menubar)
        
        fileMenu = tk.Menu(menubar)
        fileMenu.add_command(label="Add Excel File(s)", command=self.get_excel_path)
        fileMenu.add_command(label="Add Project File", command=self.get_project_path)
        fileMenu.add_command(label="Exit", command=self.master.destroy)

        menubar.add_cascade(label="File", menu=fileMenu)
        
        tk.Label(root, text="Excel Files:").grid(column=0, row=0)
        self.excel_path_entry = tk.Entry(root)
        self.excel_path_entry.grid(column=1, row=0)

        tk.Label(root, text="Project File:").grid(column=0, row=1)
        self.project_path_entry = tk.Entry(root)
        self.project_path_entry.grid(column=1, row=1)

        tk.Label(tab1, text="x-axis Format").grid(column=0, row=0)
        self.x_format = tk.StringVar(tab1)
        self.x_format.set("Absolute Time")
        tk.OptionMenu(tab1, self.x_format, *["Cycle number", 
                "Absolute Time", "Relative Time"]).grid(column=1, row=0)
        
        tk.Label(tab1, text="y-axis Format").grid(column=0, row=1)
        self.y_format = tk.StringVar(tab1)
        self.y_format.set("Concentration")
        tk.OptionMenu(tab1, self.y_format, *["Concentration", "Raw signal intensities"]).grid(column=1, row=1)

        tk.Label(tab2, text="y-axis Format").grid(column=0, row=1)
        self.y_format = tk.StringVar(tab1)
        self.y_format.set("Concentration")
        tk.OptionMenu(tab2, self.y_format, *["Concentration", "Raw signal intensities"]).grid(column=1, row=1)

        

    def get_excel_path(self):
        p = tk.filedialog.askopenfilenames(parent=root, title='Choose data files', 
        initialdir="/home/jgb509/Documents/CRM/Data")
        self.excel_path_entry.insert("0", p)
        
    def get_project_path(self):
        p = tk.filedialog.askopenfilenames(parent=root, title='Choose data files', 
        initialdir="/home/jgb509/Documents/CRM/Data")
        self.project_path_entry.insert("0", p)
    
    def load_data(self):
        print("Loading data")

    def plot_mass_scan(self):
        print("Plotting mass scan")

    def plot_time_series(self):
        print("Plotting time series")

root = tk.Tk()
app = Application(master=root)
app.mainloop()