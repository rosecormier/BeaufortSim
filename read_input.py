"""
This script reads in the data from input.txt for use by driver.py.
It saves the input parameters to a log file and prints them to the console.

R. Cormier, F. Poulin, 2023
"""

import numpy as np
#import os, sys

class Parameters():
    
    #Temporal parameters
    
    time_ave_type = None
    initial_month = None
    final_month = None
    initial_day = None 
    final_day = None
    season_start = None   
    season_end = None       

    #Spatial parameters
    
    lat_res = None      
    lon_res = None
    plot_plane_type = None
    depth_range = None
    lat_range = None  
    lon_range = None 
    # interp_type   = None #TBA
    
    #Field information
    
    # scalar: TBA
    # vector: TBA

    #Directory information
    
    logs_folder = None
    data_folder_primary = None
    data_folder_secondary = None
    visualization_folder = None
    
    clear_data_files = None #Flag to clear primary datafiles

param_data = Parameters

f = open('input.txt', 'r') #Read input file

print("Reading in parameters...")

for line in f:

    var_len = 0
    line_len = len(line)

    #Iterate over non-empty, non-commented lines
    if (line_len > 0) and (line[0] != "#"):

        for char in line:
            if char == '=':
                break
            else:
                var_len = var_len + 1

        var = line[0:var_len].strip() #Isolate variable name

        try:
            val = float(line[var_len+1:len_len-1].strip())
            val = np.array(np.matrix(val)).ravel()
        except:
            val = line[var_len+1:line_len].strip()

        setattr(param_data, var, val) #Save input data to variable

f.close

#Console output

print("Time-averaging type: ", param_data.time_ave_type)
print("Initial month: ", param_data.initial_month, "; Final month: ", param_data.final_month)
print("Initial day: ", param_data.initial_day, "; Final day: ", param_data.final_day)
print("Season start: ", param_data.season_start, "; Season end: ", param_data.season_end)

print("Resolution (latitude): ", param_data.lat_res, "; Resolution (longitude): ", param_data.lon_res)
print("Plane type: ", param_data.plot_plane_type)
print("Depth range: ", param_data.depth_range)
print("Latitude range: ", param_data.lat_range)
print("Longitude range: ", param_data.lon_range)
#print("Interpolation type: ", param_data.interp_type)

#Field info - TBA

print("Log-file directory: ", param_data.logs_folder)
print("Primary-datafile directory: ", param_data.data_folder_primary)
print("Secondary-datafile directory: ", param_data.data_folder_secondary)
print("Visualization directory: ", param_data.visualization_folder)

print("Clear primary datafiles after use: ", param_data.clear_data_files)

# save log file