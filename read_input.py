"""
This script reads in the data from input.txt for use by driver.py.
[TBA]

R. Cormier, F. Poulin, 2023
"""

import os
import sys

from os.path import expanduser, join

##############################

class Parameters():
    
    def __init__(self):
    
        #Ocean properties
        
        rho_ref = None
        nu_E = None 

        #Temporal parameters

        time_ave_type = None
        initial_month = None
        initial_year = None
        final_month = None
        final_year = None
        initial_day = None 
        final_day = None
        season_start = None   
        season_end = None
        time_kwargs = None

        #Spatial parameters

        plot_plane_type = None
        lat_res = None 
        lon_res = None
        depth_index_range = None
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
    
    def log_and_update_types(self, experiment_number):
    
        """
        Given an instance of Parameters from an input file, updates attribute 
        types to be compatible with driver.py.
        Records the input parameters in a logfile.
        """

        #Change directories from strings to paths

        self.logs_folder = join(".", self.logs_folder)
        self.visualization_folder = join(".", self.visualization_folder)
        homedir = expanduser('~')
        sys.path.append(join(homedir, 'ECCOv4-py'))
        self.data_folder_primary = join(homedir, self.data_folder_primary, 
                                        "ECCO_V4r4_PODAAC")
        self.data_folder_secondary = join(".", self.data_folder_secondary)
        
        for directory in [self.logs_folder, self.data_folder_secondary, self.visualization_folder]:
            if not os.path.exists(directory):
                os.makedirs(directory) #Make if nonexistent
                
        #Create the log file
        logfile = join(self.logs_folder, 
                       "logfile_{}.txt".format(str(experiment_number))) 
        
        f = open(logfile, "w")
        #When ready to use the script in full, switch "w" to "x"
        
        #Log input ocean properties
        f.write("rho_ref (kg/m^3) = " + self.rho_ref + "\n" + "nu_E (m^2/s) = " + 
                self.nu_E + "\n\n")
        
        #Log temporal inputs
        
        f.write("time_ave_type = " + self.time_ave_type + "\n")
        f.write("initial_month = " + self.initial_month + "\n" 
                + "initial_year = " + self.initial_year + "\n")
        f.write("final_month = " + self.final_month + "\n" + "final_year = " 
                + self.final_year + "\n\n")
        
        if self.time_ave_type == "daily":
            f.write("initial_day = " + self.initial_day + "\n" + "final_day = "
                    + self.final_day + "\n\n")
            self.time_kwargs = [self.initial_day, self.final_day]
        elif self.time_ave_type == "seasonal":
            f.write("season_start = " + self.season_start + "\n" 
                    + "season_end = " + self.season_end + "\n\n")
            self.time_kwargs = [self.season_start, self.season_end]
            
        #Log spatial inputs

        f.write("plot_plane_type = " + self.plot_plane_type + "\n\n")

        if self.plot_plane_type == "depth_index_const":
            f.write("depth = " + self.depth_index_range[0] + "\n")
            f.write("lat_res = " + self.lat_res + "\n" + "lon_res" 
                    + self.lon_res + "\n")
            f.write("latmin = " + self.lat_range[0] + "\n" + "latmax = " 
                    + self.lat_range[1] + "\n")
            f.write("lonmin = " + self.lon_range[0] + "\n" + "lonmax = " 
                    + self.lon_range[1] + "\n")
    
        elif self.plot_plane_type == "latitude_const":
            f.write("lat = " + self.lat_range[0] + "\n")
            f.write("lon_res = " + self.lon_res + "\n")
            f.write("depthmin = " + self.depth_index_range[0] + "\n" 
                    + "depthmax = " + self.depth_index_range[1] + "\n")
            f.write("lonmin = " + self.lon_range[0] + "\n" + "lonmax = " 
                    + self.lon_range[1] + "\n")
    
        elif plot_plane_type == "longitude_const":
            f.write("lon = " + self.lon_range[0] + "\n")
            f.write("lat_res = " + self.lat_res + "\n")
            f.write("depthmin = " + self.depth_index_range[0] + "\n" 
                    + "depthmax = " + self.depth_index_range[1] + "\n")
            f.write("latmin = " + self.lat_range[0] + "\n" + "latmax = " 
                    + self.lat_range[1] + "\n")
        
        #Log field inputs

        f.write("scalar_fields = ")
        for i in range(len(self.scalar_fields)-1):
            f.write(self.scalar_fields[i] + ", ")
        f.write(self.scalar_fields[-1] + "\n")
        
        f.write("vector_fields = ")
        if self.vector_fields == ['']:
            self.vector_fields = []
        else:
            for i in range(len(self.vector_fields)-1):
                f.write(self.vector_fields[i] + ", ")
        f.write(self.vector_fields[-1])
        
        return logfile
    
##############################

def parse_input_file(file, param_data):
    
    """
    Given an opened file and an instance of Parameters, parses file contents 
    and updates attributes of Parameters.
    """
    
    for line in file:

        var_len = 0
        line_len = len(line)

        #Iterate over non-empty, non-commented lines
        if (line_len > 1) and (line[0] != "#"):

            for char in line: #Isolate variable name from rest of line
                if char == '=':
                    break
                else:
                    var_len = var_len + 1

            #Strip leading/trailing whitespace from variable name
            var = line[0:var_len].strip()
            #Strip leading/trailing whitespace from remainder of line
            line = line[var_len+1:].strip()

            val_len = 0

            for char in line: 
                if char == '#': 
                    break #Isolate parameter value from trailing comment, if any
                else:
                    val_len = val_len + 1

            #Strip leading/trailing whitespace from parameter value
            val = line[0:val_len].strip()
            #Split into a list of strings
            val_list = val.split(", ") 

            if (len(val_list) > 1 or var == "scalar_fields" 
                or var == "vector_fields"):
                #If list of values (or in the cases of scalar/vector fields), 
                #save list (as strings) to variable
                setattr(param_data, var, val_list) 
            
            elif (len(val_list) == 1 and var != "scalar_fields" 
                and var != "vector_fields"): 
                #If just one value, save value (as string) to variable
                setattr(param_data, var, val) 
                
    return param_data

##############################
    
def main():
    
    param_data = Parameters()

    f = open('input.txt', 'r') #Open input file in read-only mode
    param_data = parse_input_file(f, param_data) #Parse the input file
    f.close
    
    experiment_number = 0 
    #Hardcoded for now; eventually will need to automatically produce a unique 
    #number for each run
    
    #Update parameter types (returns path to logfile)
    logfile = param_data.log_and_update_types(experiment_number)
    
    return logfile, param_data
    
##############################

if __name__ == "__main__":
    main()