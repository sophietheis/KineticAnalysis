import os
import time
import threading

import numpy as np
import pandas as pd

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from kinetic_function import (generate_profile,
                              generate_track,
                              single_track_analysis,
                              )


class ParameterApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Kinetic Analysis")
        self.root.geometry("1500x1000")

        # Bind the window close button to a custom close method
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)

        # Create Notebook (Tabs)
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(expand=True, fill="both")

        # Tab 1: Input
        self.tab1 = ttk.Frame(self.notebook)
        self.notebook.add(self.tab1, text="Generate tracks")

        # Tab 2: Result
        self.tab2 = ttk.Frame(self.notebook)
        self.notebook.add(self.tab2, text="Analyse tracks")

        # Elements for Tab 1
        self.create_input_tab()

        # Elements for Tab 2
        self.create_result_tab()

        # Initialise a list to hold CSV file names
        self.csv_files = []

    def create_input_tab(self):
        # Button to choose directory
        self.directory = None
        self.dir_button = tk.Button(self.tab1, text="Choose Directory to store output",
                                    command=lambda: self.choose_directory("tab1"))

        print(self.directory)
        self.dir_button.grid(row=0, column=0, columnspan=2, pady=10)

        # Create labels and entry widgets for 10 parameters using grid layout
        self.parameters_names = ["Protein length (aa)",
                                 "Suntag length (aa)",
                                 "Number of suntag",
                                 "Fluorescence one suntag",
                                 "Translation rate",
                                 "Binding rate",
                                 "Retention time",
                                 "Suntag position (begin or end)",
                                 "Number of tracks",
                                 "File name to save"]
        default_values = [490, 796, 32, 4, 24, 0.05, 0, "begin", 100, "datas"]
        self.parameters_input = []
        for i in range(len(self.parameters_names)):
            label = tk.Label(self.tab1, text=self.parameters_names[i])
            label.grid(row=i+1, column=0, padx=10, pady=5, sticky=tk.W)
            entry = tk.Entry(self.tab1)
            entry.insert(0, default_values[i])
            entry.grid(row=i+1, column=1, padx=10, pady=5, sticky=tk.W)
            self.parameters_input.append(entry)

        # Button to show plot
        self.plot_button = tk.Button(self.tab1, text="Show Profile", command=self.show_profile)
        self.plot_button.grid(row=11, column=0, columnspan=2, pady=10)

        # Button to start_generate_track
        self.confirm_button = tk.Button(self.tab1, text="Start generate tracks", command=self.start_generate_track)
        self.confirm_button.grid(row=12, column=0, columnspan=2, pady=10)

        # Button to save plot
        self.save_plot_button = tk.Button(self.tab1, text="Save Plot", command=self.save_plot)
        self.save_plot_button.grid(row=13, column=2, columnspan=2, pady=10)

        # Progress Bar
        self.progress = ttk.Progressbar(self.tab1, orient="horizontal", length=200, mode="determinate")
        self.progress.grid(row=14, column=0, columnspan=2, pady=10)

        # Frame for displaying the plot
        self.plot_frame = tk.Frame(self.tab1)
        self.plot_frame.grid(row=0, column=2, rowspan=13, padx=10, pady=5, sticky=tk.N)

    def create_result_tab(self):
        # Button to choose directory
        self.directory_result = None
        self.dir_button_result = tk.Button(self.tab2, text="Choose directory where your file is",
                                           command=lambda:self.choose_directory("tab2"))
        self.dir_button_result.grid(row=0, column=0, columnspan=2, pady=10)

        # Listbox to display CSV files
        self.file_listbox = tk.Listbox(self.tab2, width=50, height=10)
        self.file_listbox.grid(row=1, column=0, columnspan=2, pady=10)

        # Button to select a file from the Listbox
        self.select_button = tk.Button(self.tab2, text="Select File", command=self.select_file)
        self.select_button.grid(row=2, column=0, columnspan=2, pady=10)

        # Create labels and entry widgets for 10 parameters using grid layout
        self.parameters_names_analysis = ["dt",
                                          "protein length (aa)",
                                          "File name to save"]
        default_values = [3, 800, "datas_results"]
        self.parameters_analysis = []
        for i in range(len(self.parameters_names_analysis)):
            label = tk.Label(self.tab2, text=self.parameters_names_analysis[i])
            label.grid(row=i + 3, column=0, padx=10, pady=5, sticky=tk.W)
            entry = tk.Entry(self.tab2)
            entry.insert(0, default_values[i])
            entry.grid(row=i + 3, column=1, padx=10, pady=5, sticky=tk.W)
            self.parameters_analysis.append(entry)

        # Button to start analysis track
        self.confirm_analyse_button = tk.Button(self.tab2, text="Start analyse tracks", command=self.start_analyse_track)
        self.confirm_analyse_button.grid(row=10, column=0, columnspan=2, pady=10)

        # Progress Bar
        self.progress_analyse = ttk.Progressbar(self.tab2, orient="horizontal", length=200, mode="determinate")
        self.progress_analyse.grid(row=12, column=0, columnspan=2, pady=10)

    def choose_directory(self, tab):
        directory = filedialog.askdirectory()
        if tab == "tab1":
            self.directory = directory
        elif tab == "tab2":
            self.directory_result = directory
            self.load_csv_files()
        if directory and tab !="tab2":
            messagebox.showinfo("Directory Selected", f"Directory chosen: {directory}")

    def load_csv_files(self):
        """Load CSV files from the selected directory into the Listbox."""
        # Clear the Listbox
        self.file_listbox.delete(0, tk.END)
        self.csv_files = []  # Clear previous file list

        # List all CSV files in the directory
        if self.directory_result:  # Ensure a directory is selected
            for file in os.listdir(self.directory_result):
                if file.endswith('.csv'):
                    self.csv_files.append(file)
                    self.file_listbox.insert(tk.END, file)  # Add file to the Listbox

    def select_file(self):
        """Select a file from the Listbox and display a message."""
        selected_indices = self.file_listbox.curselection()
        if selected_indices:
            self.selected_file_analyse = self.file_listbox.get(selected_indices[0])
            # messagebox.showinfo("File Selected", f"You selected: {self.selected_file_analyse}")
        else:
            messagebox.showwarning("Warning", "No file selected.")

    def show_profile(self):
        # Clear any existing plot
        for widget in self.plot_frame.winfo_children():
            widget.destroy()

        # Retrieve parameter values
        param_values = [entry.get() for entry in self.parameters_input]

        # Check if all parameters are filled
        if all(param_values):
            # Create graphes
            self.fig = Figure(figsize=(6, 5), dpi=100)
            ax = self.fig.add_subplot(211)
            x, y = generate_profile(float(param_values[0]),
                                    float(param_values[1]),
                                    float(param_values[2]),
                                    float(param_values[3]),
                                    float(param_values[4]),
                                    float(param_values[5]),
                                    float(param_values[6]),
                                    param_values[7],
                                    )
            ax.plot(x, y)
            ax.set_title("One protein fluorescence profile")
            ax.set_xlabel("Time")
            ax.set_ylabel("Fluorescence")

            ax = self.fig.add_subplot(212)
            if param_values[7] == "begin":
                suntag_pos = 0
            else:
                suntag_pos = -1
            x, y, _ = generate_track(float(param_values[0]),
                                     float(param_values[1]),
                                     float(param_values[2]),
                                     float(param_values[3]),
                                     float(param_values[4]),
                                     float(param_values[5]),
                                     float(param_values[6]),
                                     suntag_pos)
            ax.plot(x, y)
            ax.set_title("Fluorescence profile")
            ax.set_xlabel("Time")
            ax.set_ylabel("Fluorescence")
            # Embed the plot in the Tkinter window
            canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
            canvas.draw()
            canvas.get_tk_widget().pack()

        else:
            messagebox.showwarning("Warning", "Please fill all parameters.")

    def start_generate_track(self):
        # Start progress bar in a separate thread
        threading.Thread(target=self.run_progress_and_start_generate_track).start()

    def run_progress_and_start_generate_track(self):
        # Retrieve parameter values
        param_values = [entry.get() for entry in self.parameters_input]

        # Check if all parameters are filled
        if all(param_values) and self.directory:
            # Initialize progress bar
            self.progress["maximum"] = int(param_values[8])
            self.progress["value"] = 0

            # Generate all tracks and save it
            first_time = True
            if param_values[7] == "begin":
                suntag_pos = 0
            else:
                suntag_pos = -1
            for i in range(int(param_values[8])):
                x_global, y_global, y_start_prot = generate_track(float(param_values[0]),
                                                                  float(param_values[1]),
                                                                  float(param_values[2]),
                                                                  float(param_values[3]),
                                                                  float(param_values[4]),
                                                                  float(param_values[5]),
                                                                  float(param_values[6]),
                                                                  suntag_pos)
                if first_time:
                    datas = pd.DataFrame({"FRAME": x_global,
                                          "MEAN_INTENSITY_CH1": y_global,
                                          "TRACK_ID": i,
                                          "RETENTION_TIME": float(param_values[6]),
                                          })
                    first_time = False
                else:
                    datas = pd.concat([datas,
                                       pd.DataFrame({"FRAME": x_global,
                                                     "MEAN_INTENSITY_CH1": y_global,
                                                     "TRACK_ID": i,
                                                     "RETENTION_TIME": float(param_values[6]),
                                                     })], ignore_index=True)
                self.progress["value"] = i  # Update the progress bar
                self.root.update_idletasks()

            datas.to_csv(os.path.join(self.directory, param_values[9] + ".csv"))

            self.progress["value"] = i + 1  # Update the progress bar
            self.root.update_idletasks()

        else:
            messagebox.showwarning("Warning", "Please fill all parameters and select a directory.")

    def start_analyse_track(self):
        # Start progress bar in a separate thread
        threading.Thread(target=self.run_progress_and_start_analyse_track).start()

    def run_progress_and_start_analyse_track(self):
        # Retrieve parameter values
        param_values = [entry.get() for entry in self.parameters_analysis]

        # Check if all parameters are filled
        if all(param_values) and self.directory_result:

            # Read csv file
            datas = pd.read_csv(os.path.join(self.directory_result, self.selected_file_analyse),
                                index_col="Unnamed: 0")
            dt = float(param_values[0])
            t = dt/0.1
            prot_length = float(param_values[1])
            nb_track = len(np.unique(datas["TRACK_ID"]))

            # Initialize progress bar
            self.progress_analyse["maximum"] = nb_track
            self.progress_analyse["value"] = 0

            first_time = True
            # Analyse all tracks and save it
            for i in range(nb_track):
                datas2 = datas[(datas["TRACK_ID"] == i)][::int(t)]

                (x,
                 y,
                 x_auto,
                 y_auto,
                 elongation_r,
                 translation_init_r,
                 perr) = single_track_analysis(datas2,
                                               i,
                                               delta_t=dt,
                                               protein_size=prot_length,
                                               normalise_intensity=1,
                                               normalize_auto=True,
                                               mm=None,
                                               lowpass_=False,
                                               cutoff=100,
                                               rtol=1e-1,
                                               method="linear",
                                               force_analysis=True,
                                               first_dot=True,
                                               simulation=True)
                if first_time:
                    results = pd.DataFrame({"elongation_r": elongation_r,
                                            "init_translation_r": translation_init_r,
                                            "dt": dt,
                                            "id": i,},
                                           index=[0])
                    first_time = False

                else:
                    results = pd.concat([results,
                                         pd.DataFrame({"elongation_r": elongation_r,
                                                       "init_translation_r": translation_init_r,
                                                       "dt": dt,
                                                       "id": i,}, index=[0])
                                         ], ignore_index=True)

                self.progress_analyse["value"] = i  # Update the progress bar
                self.root.update_idletasks()
            results.to_csv(os.path.join(self.directory_result, param_values[2]+".csv"))

            self.progress_analyse["value"] = i + 1  # Update the progress bar
            self.root.update_idletasks()

        else:
            messagebox.showwarning("Warning", "Please fill all parameters and select a directory.")

    def save_plot(self):
        """Save the currently displayed plot to a file."""
        if self.fig is not None:
            # Ask the user for a file location to save the plot
            file_path = filedialog.asksaveasfilename(
                defaultextension=".png",
                filetypes=[("PNG files", "*.png"), ("EPS files", "*.eps"), ("All Files", "*.*")],
                title="Save Plot As"
            )
            if file_path:
                # Save the plot to the specified file path
                try:
                    self.fig.savefig(file_path)
                    messagebox.showinfo("Success", f"Plot saved as {file_path}")
                except Exception as e:
                    messagebox.showerror("Error", f"Failed to save plot: {str(e)}")
        else:
            messagebox.showwarning("Warning", "No plot to save. Please generate a plot first.")

    def on_close(self):
        """Handle the window close event."""
        # Confirm before closing the window
        if messagebox.askokcancel("Quit", "Do you really want to quit?"):
            self.root.destroy()  # Properly close the application


if __name__ == "__main__":
    root = tk.Tk()
    app = ParameterApp(root)
    root.mainloop()
