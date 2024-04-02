import nmrglue as ng
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from tkinter import *
import os
import zipfile
import shutil
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from scipy.signal import find_peaks
from pandastable import Table, TableModel
from scipy.integrate import simps
import sv_ttk

# Global variables
selected_pdata_dirs = []
tsp_concentration = 0
binning_step = 0.05
spectra = []
peak_limits = pd.DataFrame()

def find_pdata_directories(root_dir):
    pdata_dirs = []
    for dirpath, _, filenames in os.walk(root_dir):
        if "/pdata/1" in dirpath:
            pdata_dirs.append(dirpath)
    return pdata_dirs

def browse_directory():
    global selected_pdata_dirs
    global spectra
    global peak_limits
    file_or_dir = messagebox.askquestion("Zipped File or Directory", "Are you selecting a zipped file?", icon='question')
    if file_or_dir == 'yes':
            root_dir = filedialog.askopenfilename(title="Select the zipped file containing processed Bruker NMR data")
            if not root_dir: # User closed the dialog without selecting a file
                messagebox.showerror("Error", "No zipped file selected. Please select a valid file.")
                return
            if not root_dir.endswith('.zip'):
                messagebox.showerror("Error", "Selected file is not a zip file. Please select a valid zip file.")
                return
            with zipfile.ZipFile(root_dir, 'r') as zip_ref:
                zip_ref.extractall('Zipped file')
                root_dir = 'Zipped file'
    else:
        root_dir = filedialog.askdirectory(title="Select the root directory containing processed Bruker NMR data")
        if not root_dir: # User closed the dialog without selecting a directory
            messagebox.showerror("Error", "No directory selected. Please select a valid directory.")
            return
    selected_pdata_dirs = [root_dir]

    peak_limits = pd.read_excel("peak_limits.xlsx")

    spectra = []
    for root_dir in selected_pdata_dirs:
        pdata_dirs = find_pdata_directories(root_dir)
        for pdata_dir in pdata_dirs:
            try:
                dic, data = ng.bruker.read_pdata(pdata_dir, scale_data=True)
                udic = ng.bruker.guess_udic(dic, data)
                uc = ng.fileiobase.uc_from_udic(udic)
                ppm_scale = uc.ppm_scale()
                spectra.append((ppm_scale, data))
            except OSError:
                continue
    selected_dirs_label.config(text=f"Path Chosen: {root_dir}")

    plot_spectra(spectra, peak_limits)

def process_selected_dirs_concentration():
    global selected_pdata_dirs, tsp_concentration

    if not selected_pdata_dirs:
        messagebox.showerror("Error", "Please select a data directory.")
        return

    try:
        tsp_concentration = float(concentration_entry.get())
    except ValueError:
        messagebox.showerror("Error", "Invalid concentration value. Please enter a valid number.")
        return

    results_concentration = []
    results_area = []
    peak_identities = set()

    for root_dir in selected_pdata_dirs:
        pdata_dirs = find_pdata_directories(root_dir)
        if not pdata_dirs:
            messagebox.showwarning("No pdata/1 Directories", f"No pdata/1 directories found in '{root_dir}'. Skipping.")
            continue

        for pdata_dir in pdata_dirs:
            try:
                dic, data = ng.bruker.read_pdata(pdata_dir, scale_data=True)
            except OSError:
                continue

            udic = ng.bruker.guess_udic(dic, data)
            uc = ng.fileiobase.uc_from_udic(udic)
            ppm_scale = uc.ppm_scale()

            peak_limits = pd.read_excel("peak_limits.xlsx")

            ref_name = peak_limits.at[0, 'Peak identity']
            ref_start = peak_limits.at[0, 'ppm start']
            ref_end = peak_limits.at[0, 'ppm end']
            ref_num_protons = peak_limits.at[0, '# protons']

            ref_min_index = np.abs(ppm_scale - ref_start).argmin()
            ref_max_index = np.abs(ppm_scale - ref_end).argmin()
            if ref_min_index > ref_max_index:
                ref_min_index, ref_max_index = ref_max_index, ref_min_index

            ref_peak = data[ref_min_index:ref_max_index + 1]
            ref_area = ref_peak.sum()

            for index, row in peak_limits.iloc[1:].iterrows():
                name = row['Peak identity']
                peak_identities.add(name)
                start = row['ppm start']
                end = row['ppm end']
                num_protons_peak = row['# protons']

                min_index = np.abs(ppm_scale - start).argmin()
                max_index = np.abs(ppm_scale - end).argmin()
                if min_index > max_index:
                    min_index, max_index = max_index, min_index

                peak = data[min_index:max_index + 1]
                peak_area = peak.sum()
                results_area.append({'Name': name, 'Area': peak_area, 'Parent File Path': pdata_dir})

                if tsp_concentration != 0:
                    concentration = (peak_area / ref_area) * tsp_concentration * ref_num_protons / num_protons_peak
                    results_concentration.append({'Name': name, 'Concentration': concentration, 'Parent File Path': pdata_dir})

    if not results_concentration:
        messagebox.showerror("No Data", "No pdata/1 directories were processed.")
        return

    concentrations_df = pd.DataFrame(results_concentration)
    concentrations_df = concentrations_df.pivot(index='Parent File Path', columns='Name', values='Concentration')

    areas_df = pd.DataFrame(results_area)
    areas_df = areas_df.pivot(index='Parent File Path', columns='Name', values='Area')

    with pd.ExcelWriter('nmr_analysis_results.xlsx') as writer:
        concentrations_df.to_excel(writer, sheet_name='Absolute Concentrations')
        areas_df.to_excel(writer, sheet_name='Areas')

    messagebox.showinfo("Processing Complete", "Check the results in 'nmr_analysis_results.xlsx'")

    # Clean up temporary directory from unzipping
    for dir in selected_pdata_dirs:
        if 'temp' in dir:
            shutil.rmtree(dir)

def process_selected_dirs_binning():
    global selected_pdata_dirs

    if not selected_pdata_dirs:
        messagebox.showerror("Error", "Please select a data directory.")
        return

    try:
        binning_step_size = float(binning_entry.get())
    except ValueError:
        messagebox.showerror("Error", "Invalid binning step size. Please enter a valid number.")
        return

    all_bins = []
    all_ppm_scales = []

    for root_dir in selected_pdata_dirs:
        pdata_dirs = find_pdata_directories(root_dir)
        if not pdata_dirs:
            messagebox.showwarning("No pdata/1 Directories", f"No pdata/1 directories found in '{root_dir}'. Skipping.")
            continue

        for pdata_dir in pdata_dirs:
            try:
                dic, data = ng.bruker.read_pdata(pdata_dir, scale_data=True)
            except OSError:
                continue

            udic = ng.bruker.guess_udic(dic, data)
            uc = ng.fileiobase.uc_from_udic(udic)
            ppm_scale = uc.ppm_scale()
            all_ppm_scales.append(ppm_scale)
            
    min_ppm = min([ppm.min() for ppm in all_ppm_scales])
    max_ppm = max([ppm.max() for ppm in all_ppm_scales])

    bin_edges = np.arange(min_ppm, max_ppm, binning_step_size)

    for ppm_scale in all_ppm_scales:
        bin_indices = np.digitize(ppm_scale, bin_edges)
        binned_spectrum = [data[bin_indices == i].sum() for i in range(1, len(bin_edges))]

        all_bins.append(binned_spectrum)

    binning_df = pd.DataFrame(all_bins, columns=bin_edges[:-1], index=[os.path.dirname(pdata_dir) for pdata_dir in pdata_dirs])

    binning_df = binning_df.T

    with pd.ExcelWriter('nmr_binning_results.xlsx') as writer:
        binning_df.to_excel(writer, sheet_name='Binning')

    messagebox.showinfo("Binning Complete", "Check the results in 'nmr_binning_results.xlsx'")

    # Clean up temporary directories
    for dir in selected_pdata_dirs:
        if 'temp' in dir:
            shutil.rmtree(dir)

def plot_spectra(spectra, peak_data):
    fig.clear() # Clear the entire figure
    for ax in fig.axes:
        ax.cla()  # clear the plot
    ax = fig.add_subplot(111)

    max_y_values = {} 
    min_y_values = {} # Dictionary to store the min and max y value for each peak across all spectra

    for i, (ppm_scale, data) in enumerate(spectra):
        ax.plot(ppm_scale, data)
        
        for index, row in peak_data.iterrows():
            name = row['Peak identity']
            start = row['ppm start']
            end = row['ppm end']
            mid_point = (start + end) / 2  # Calculate midpoint for label placement

            # Get indices of start and end points
            start_index = np.abs(ppm_scale - start).argmin()
            end_index = np.abs(ppm_scale - end).argmin()
            if start_index > end_index:
                start_index, end_index = end_index, start_index

            # Get the maximum and minimum y value in the range of the peak (start to end)
            max_y_in_range = data[start_index:end_index + 1].max()
            min_y_in_range = data[start_index:end_index + 1].min()
            
            # If this peak's maximum y value is the largest one found so far, store it
            if name not in max_y_values or max_y_in_range > max_y_values[name]:
                max_y_values[name] = max_y_in_range
            # If this peak's minimum y value is the largest one found so far, store it
            if name not in min_y_values or min_y_in_range < min_y_values[name]:
                min_y_values[name] = min_y_in_range

            # Add shading for the integration area
            ax.fill_between(ppm_scale[start_index:end_index + 1], data[start_index:end_index + 1], color='blue', alpha=0.25)

    for index, row in peak_data.iterrows():
        name = row['Peak identity']
        start = row['ppm start']
        end = row['ppm end']
        ax.plot([start, end], [min_y_values[name], min_y_values[name]], color='black', linewidth=1)  # Draw bottom line for integration region
        ax.plot([start, end], [max_y_values[name], max_y_values[name]], color='black', linewidth=1) # Draw top line

        mid_point = (start + end) / 2  # Horizontal label placement
        offset = max_y_values[name] * 0.01 # Vertical label placement
        ax.text(mid_point, max_y_values[name] + offset, name, ha='center', va='bottom', color='black', fontsize=10, rotation=90) # Print text vertically
        
    ax.set_xlabel("Chemical Shift (ppm)", fontsize=16)
    ax.set_ylabel("Intensity", fontsize=16)
    ax.invert_xaxis()
    ax.format_coord = lambda x, y: ""
    ax.autoscale_view()
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=1)

selected_pdata_dirs = []
tsp_concentration = 0
binning_step = 0.05

###GUI###
root = tk.Tk()
root.title("Bruker NMR Analysis")
root.iconbitmap('/icon.ico')
root.geometry("1000x750")

title_main = tk.Label(root, text="Bruker NMR Analysis", font=('','24','bold'))
title_main.pack()

###TABS###
tabControl = ttk.Notebook(root)

tab_selectData = ttk.Frame(tabControl)
tab_concentration = ttk.Frame(tabControl)
tab_binning = ttk.Frame(tabControl)

tabControl.add(tab_selectData, text ='Select Dataset')
tabControl.add(tab_concentration, text ='Absolute Concentrations')
tabControl.add(tab_binning, text ='Binning')
tabControl.pack(expand=1, fill="none")

# DATA SELECTION START
select_text1 = tk.Label(tab_selectData, text="Select Processed Bruker NMR dataset: ", font=('','14','bold'))
select_text1.grid(row=0, sticky=W)

select_button = tk.Button(tab_selectData, text="Choose file", command=browse_directory)
select_button.grid(row=0, column=1, sticky=W)

selected_dirs_label = tk.Label(tab_selectData, text="")
selected_dirs_label.grid(row=1, sticky=W)

select_text2 = tk.Label(tab_selectData, text="Information from peak_limits will be overlayed")
select_text2.grid(row=2, sticky=W)

select_text3 = tk.Label(tab_selectData, text="Edit peak_limits as needed")
select_text3.grid(row=3, sticky=W)

select_text4 = tk.Label(tab_selectData, text="peak_limits.xlsx contains metabolite peak information")
select_text4.grid(row=4, sticky=W)

# CONCENTRATION START
title_2 = tk.Label(tab_concentration, text="Enter Chemical Shift Reference Concentration: ", font=('','14','bold'))
title_2.grid(row=0, sticky=W)

concentration_label = tk.Label(tab_concentration, text="Enter the Chemical Shift Reference Concentration for the Samples: ")
concentration_label.grid(row=1, sticky=E)
concentration_entry = tk.Entry(tab_concentration)
concentration_entry.grid(row=1, column=1, sticky=W)
concentration_label = tk.Label(tab_concentration, text="uM")
concentration_label.grid(row=1, column=1, sticky=E)

title_3 = tk.Label(tab_concentration, text="Process Area and Absolute Concentrations", font=('','14','bold'))
title_3.grid(row=3, sticky=W)

process_concentration_button = tk.Button(tab_concentration, text="Begin", command=process_selected_dirs_concentration)
process_concentration_button.grid(row=4, sticky=W)

# BINNING START
title_4 = tk.Label(tab_binning, text="Perform Binning and Calculate Areas: ", font=('','14','bold'))
title_4.grid(row=0, sticky=W)

binning_label = tk.Label(tab_binning, text="Enter the Step Size, Width of each Bin: ")
binning_label.grid(row=1,sticky=E)
binning_entry = tk.Entry(tab_binning)
binning_entry.grid(row=1,column=1,sticky=W)
binning_label = tk.Label(tab_binning, text="ppm")
binning_label.grid(row=1,column=1,sticky=E)

binning_button_headline = tk.Label(tab_binning, text="Process Binning Areas", font=('','14','bold'))
binning_button_headline.grid(row=2,sticky=W)
process_binning_button = tk.Button(tab_binning, text="Begin", command=process_selected_dirs_binning)
process_binning_button.grid(row=3, sticky=W)

###PLOT###
plotFrame = Frame(root)
fig = plt.Figure(figsize=(5, 5), dpi=100)
canvas = FigureCanvasTkAgg(fig, master=root)
toolbar = NavigationToolbar2Tk(canvas, plotFrame)

toolbar.update()
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
canvas.get_tk_widget().pack()
plotFrame.pack()

sv_ttk.set_theme("light")
root.mainloop()