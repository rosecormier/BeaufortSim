"""
This script reads in the data from input.txt for use by driver.py.
It prints the input parameters to the console.

R. Cormier, F. Poulin, 2023
"""

import numpy as np

##############################

class Parameters():
    
    #Ocean properties
    rho_ref = None
    nu_E = None 
    
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
    
    scalar_fields = None
    vector_fields = None

    #Directory information
    
    logs_folder = None
    data_folder_primary = None
    data_folder_secondary = None
    visualization_folder = None
    
    clear_data_files = None #Flag to clear primary datafiles

##############################
    
def main():
    
    param_data = Parameters

    f = open('input.txt', 'r') #Read input file

    print("Reading in parameters...")

    for line in f:

        var_len = 0
        line_len = len(line)

        #Iterate over non-empty, non-commented lines
        if (line_len > 1) and (line[0] != "#"):

            for char in line: #Isolate variable name from rest of line
                if char == '=':
                    break
                else:
                    var_len = var_len + 1

            var = line[0:var_len].strip() #Strip leading/trailing whitespace from variable name
            line = line[var_len+1:].strip() #Strip leading/trailing whitespace from remainder of line

            val_len = 0

            for char in line: #Isolate parameter value from trailing comment, if any
                if char == '#': 
                    break
                else:
                    val_len = val_len + 1

            val = line[0:val_len].strip() #Strip leading/trailing whitespace from parameter value
            val_list = val.split(", ") #Split into a list of strings

            if len(val_list) > 1 or var == "scalar_fields" or var == "vector_fields":
                setattr(param_data, var, val_list) #If list of values (or in the cases of scalar/vector fields), save list (as strings) to variable
            
            elif len(val_list) == 1 and var != "scalar_fields" and var != "vector_fields": 
                setattr(param_data, var, val) #If just one value, save value (as string) to variable

    f.close

    #Console output
    
    print("Reference density (kg/m^3):", param_data.rho_ref)
    print("Eddy viscosity (m^2/s):", param_data.nu_E)

    print("Time-averaging type:", param_data.time_ave_type)
    print("Initial month:", param_data.initial_month[1]+",", param_data.initial_month[0])
    print("Final month:", param_data.final_month[1]+",", param_data.final_month[0])
    print("Initial day:", param_data.initial_day)
    print("Final day:", param_data.final_day)
    print("First month of season:", param_data.season_start)
    print("Last month of season:", param_data.season_end)

    print("Resolution (latitude):", param_data.lat_res)
    print("Resolution (longitude):", param_data.lon_res)
    print("Plane type:", param_data.plot_plane_type)
    print("Depth range:", param_data.depth_range[0]+",", param_data.depth_range[1])
    print("Latitude range:", param_data.lat_range[0]+",", param_data.lat_range[1])
    print("Longitude range:", param_data.lon_range[0]+",", param_data.lon_range[1])
    #print("Interpolation type:", param_data.interp_type)

    print("Scalar fields:", param_data.scalar_fields)
    print("Vector fields:", param_data.vector_fields)

    print("Log-file directory:", param_data.logs_folder)
    print("Primary-datafile directory:", param_data.data_folder_primary)
    print("Secondary-datafile directory:", param_data.data_folder_secondary)
    print("Visualization directory:", param_data.visualization_folder)

    print("Clear primary datafiles after use:", param_data.clear_data_files)
    
    return param_data
    
##############################

if __name__ == "__main__":
    main()