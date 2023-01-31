import tkinter as tk
from tkinter import filedialog as fd


class MainGUI(tk.Tk):

    def __init__(self):
        super().__init__()

        self.wm_title("CA Parameters")
        self.title("CA Parameters")
        # self.geometry("500x500")

        self.ehooke_var = tk.StringVar()
        self.ehooke_var.set("None")

        self.root_var = tk.StringVar()
        self.root_var.set("None")

        self.base_var = tk.StringVar()
        self.base_var.set("")

        self.fluor1_var = tk.StringVar()
        self.fluor1_var.set("")

        self.fluor2_var = tk.StringVar()
        self.fluor2_var.set("")

        self.membrane_var = tk.StringVar()
        self.membrane_var.set("")

        self.basetype_var = tk.StringVar()
        self.basetype_var.set("Phase")

        self.pars = None

        self.init_files()
        self.init_names()
        self.init_basetype()

        tk.Button(master=self, text="ENTER", command=self.return_parameters).pack(fill='both', expand=True, side='top')

    def init_files(self):
        files_frame = tk.Frame(self)
        files_frame.pack(fill='both', expand=True, side='top')

        files_frame.columnconfigure(1, weight=1)

        ehooke_button = tk.Button(master=files_frame, text="eHooke Folder", command=self.getehookepath)
        ehooke_button.grid(column=0, row=0, sticky=tk.W + tk.E)

        ehooke_label = tk.Label(master=files_frame, textvariable=self.ehooke_var)
        ehooke_label.grid(column=1, row=0, sticky=tk.W + tk.E)

        root_button = tk.Button(master=files_frame, text="Root Folder", command=self.getrootpath)
        root_button.grid(column=0, row=1, sticky=tk.W + tk.E)

        root_label = tk.Label(master=files_frame, textvariable=self.root_var)
        root_label.grid(column=1, row=1, sticky=tk.W + tk.E)

    def init_names(self):
        names_frame = tk.Frame(self)
        names_frame.pack(fill='both', expand=True, side='top')

        fluor1_label = tk.Label(master=names_frame, text="Fluorescence file name", anchor='w')
        fluor1_label.grid(column=0, row=0, sticky=tk.W + tk.E)
        fluor1_entry = tk.Entry(master=names_frame, textvariable=self.fluor1_var)
        fluor1_entry.grid(column=1, row=0, sticky=tk.W + tk.E)

        fluor2_label = tk.Label(master=names_frame, text="Fluorescence file name (optional)", anchor='w')
        fluor2_label.grid(column=0, row=1, sticky=tk.W + tk.E)
        fluor2_entry = tk.Entry(master=names_frame, textvariable=self.fluor2_var)
        fluor2_entry.grid(column=1, row=1, sticky=tk.W + tk.E)

        memb_label = tk.Label(master=names_frame, text="Membrane file name (optional)", anchor='w')
        memb_label.grid(column=0, row=2, sticky=tk.W + tk.E)
        memb_entry = tk.Entry(master=names_frame, textvariable=self.membrane_var)
        memb_entry.grid(column=1, row=2, sticky=tk.W + tk.E)

        base_label = tk.Label(master=names_frame, text="Base file name", anchor='w')
        base_label.grid(column=0, row=3, sticky=tk.W + tk.E)
        base_entry = tk.Entry(master=names_frame, textvariable=self.base_var)
        base_entry.grid(column=1, row=3, sticky=tk.W + tk.E)

    def init_basetype(self):
        basetype_frame = tk.Frame(self)
        basetype_frame.pack(fill='both', expand=True, side='top')
        basetype_option = {"Phase": "Phase", "BF": "BF", "Membrane": "Membrane"}

        basetype_label = tk.Label(master=basetype_frame, text="Choose type of base image used")
        basetype_label.pack()

        for txt, val in basetype_option.items():
            tk.Radiobutton(master=basetype_frame, text=txt, variable=self.basetype_var, value=val) \
                .pack(fill='both', expand=True, side='left')

    def getrootpath(self):
        p = fd.askdirectory(initialdir="C:", title="Select Folder")
        self.root_var.set(p)

    def getehookepath(self):
        p = fd.askdirectory(initialdir="C:", title="Select Folder")
        self.ehooke_var.set(p)

    def return_parameters(self):
        self.pars = {"ehookefolder": self.ehooke_var.get(),
                     "rootfolder": self.root_var.get(),
                     "fluor1": self.fluor1_var.get(),
                     "fluor2": self.fluor2_var.get(),
                     "base": self.base_var.get(),
                     "basetype": self.basetype_var.get()}

        self.destroy()

if __name__ == '__main__':
    MainGUI().mainloop()
