"""
This script reads in the data from input.txt for use by driver.py.
It saves the input parameters to a log file and prints them to the console.

R. Cormier, F. Poulin, 2023
"""

f = open('input.txt', 'r')

import numpy as np
import os, sys

class Parameters():

    initial_date = None
    final_date   = None

param_data = Parameters

for line in f:

    var_len = 0
    line_len = len(line)

    if (line_len > 0) and (line[0] != "#"):

        for char in line:
            if char == '=':
                break
            else:
                var_len = var_len + 1

        var = line[0:var_len].strip()

        try:
            val = float(line[var_len+1:len_len-1].strip())
            val = np.array(np.matrix(val)).ravel()
        except:
            val = line[var_len+1:line_len].strip()

        setattr(param_data, var, val)

f.close

print(param_data.initial_date)
print(param_data.final_date)

print(param_data.time_ave_type)
print(param_data.monthly_range)

print(param_data.lat_range)
print(param_data.lon_range)
print(param_data.depth_range)
print(param_data.lat_res)
print(param_data.lon_res)
print(param_data.interp_type)

print(param_data.plot_plane_type)
print(param_data.k_plot)
print(param_data.lat_plot)
print(param_data.lon_plot)

print(param_data.plot_fields)
print(param_data.data_folder_primary)
print(param_data.data_folder_secondary)
print(param_data.visualization_folder)

# save log file