"""
This script uses the data provided in the input file (input.txt) to do all of the following:
    -Save a txt file listing the input parameters to be used;
    -Download data corresponding to specified fields;
    -Computes and saves secondary fields (skips if unnecessary or done);
    -Plots data for all fields specified and saves plots; and
    -Removes primary data files.

R. Cormier, F. Poulin, 2023
"""

import os
import sys

from os.path import expanduser, join

import read_input
import download_data

from functions_comp_data_meta import create_comp_data_file, load_comp_data_file
from NEWfunctions_ecco_general import load_ECCO_data_file
from functions_field_variables import field_is_primary

##############################

#READ (AND LOG, AS APPROPRIATE) INPUT PARAMETERS

param_data = read_input.main()

#Input directory information

logs_folder = param_data.logs_folder
data_folder_primary = param_data.data_folder_primary
data_folder_secondary = param_data.data_folder_secondary
visualization_folder = param_data.visualization_folder

#Set up directories

homedir = expanduser('~')
sys.path.append(join(homedir, 'ECCOv4-py'))
datdir_primary = join(homedir, data_folder_primary, 'ECCO_V4r4_PODAAC')
datdir_secondary = join(".", data_folder_secondary)
visdir = join(".", visualization_folder)
logdir = join(".", logs_folder)

experiment_number = 0 #Hardcoded for now; eventually will need to automatically produce a unique number for each run

if not os.path.exists(logdir):
    os.makedirs(logdir)

logfile = join(logdir, "logfile_{}.txt".format(str(experiment_number))) #Create the log file

f = open(logfile, "w") #When ready to use the script in full, switch "w" to "x"

#Ocean properties

rho_ref = param_data.rho_ref
nu_E = param_data.nu_E
f.write("rho_ref (kg/m^3) = " + rho_ref + "\n" + "nu_E (m^2/s) = " + nu_E + "\n\n")

#Temporal inputs

time_ave_type = param_data.time_ave_type
initial_year, initial_month = param_data.initial_month[0], param_data.initial_month[1]
final_year, final_month = param_data.final_month[0], param_data.final_month[1]
f.write("time_ave_type = " + time_ave_type + "\n")
f.write("initial_month = " + initial_month + "\n" + "initial_year = " + initial_year + "\n")
f.write("final_month = " + final_month + "\n" + "final_year = " + final_year + "\n\n")

if time_ave_type == "daily":
    initial_day, final_day = param_data.initial_day, param_data.final_day
    f.write("initial_day = " + initial_day + "\n" + "final_day = " + final_day + "\n\n")
elif time_ave_type == "season":
    season_start, season_end = param_data.season_start, param_data.season_end
    f.write("season_start = " + season_start + "\n" + "season_end = " + season_end + "\n\n")

#Spatial inputs
    
plot_plane_type = param_data.plot_plane_type
f.write("plot_plane_type = " + plot_plane_type + "\n\n")

if plot_plane_type == "depth_const":
    depth = param_data.depth_range[0]
    lat_res, lon_res = param_data.lat_res, param_data.lon_res
    latmin, latmax = param_data.lat_range[0], param_data.lat_range[1]
    lonmin, lonmax = param_data.lon_range[0], param_data.lon_range[1]
    f.write("depth = " + depth + "\n")
    f.write("lat_res = " + lat_res + "\n" + "lon_res" + lon_res + "\n")
    f.write("latmin = " + latmin + "\n" + "latmax = " + latmax + "\n")
    f.write("lonmin = " + lonmin + "\n" + "lonmax = " + lonmax + "\n")
elif plot_plane_type == "latitude_const":
    lat = param_data.lat_range[0]
    lon_res = param_data.lon_res
    depthmin, depthmax = param_data.depth_range[0], param_data.depth_range[1]
    lonmin, lonmax = param_data.lon_range[0], param_data.lon_range[1]
    f.write("lat = " + lat + "\n")
    f.write("lon_res = " + lon_res + "\n")
    f.write("depthmin = " + depthmin + "\n" + "depthmax = " + depthmax + "\n")
    f.write("lonmin = " + lonmin + "\n" + "lonmax = " + lonmax + "\n")
elif plot_plane_type == "longitude_const":
    lon = param_data.lon_range[0]
    lat_res = param_data.lat_res
    depthmin, depthmax = param_data.depth_range[0], param_data.depth_range[1]
    latmin, latmax = param_data.lat_range[0], param_data.lat_range[1]
    f.write("lon = " + lon + "\n")
    f.write("lat_res = " + lat_res + "\n")
    f.write("depthmin = " + depthmin + "\n" + "depthmax = " + depthmax + "\n")
    f.write("latmin = " + latmin + "\n" + "latmax = " + latmax + "\n")
    
# interp_type = #TBA

#Field inputs

scalar_fields = param_data.scalar_fields
f.write("scalar_fields = ")
for i in range(len(scalar_fields)-1):
    f.write(scalar_fields[i] + ", ")
f.write(scalar_fields[-1] + "\n")
        
vector_fields = param_data.vector_fields
f.write("vector_fields = ")
for i in range(len(vector_fields)-1):
    f.write(vector_fields[i] + ", ")
f.write(vector_fields[-1])

f.close()
    
#Flag to clear primary datafiles after use (no need to log)
clear_data_files = param_data.clear_data_files 

##############################

#DISTINGUISH PRIMARY (ECCO) FIELDS FROM SECONDARY (COMPUTED) FIELDS

primary_scalar_fields, secondary_scalar_fields = [], []

for scalar_field_name in scalar_fields:
    if field_is_primary(scalar_field_name):
        primary_scalar_fields.append(scalar_field_name)
    elif not field_is_primary(scalar_field_name):
        secondary_scalar_fields.append(scalar_field_name)
        
primary_vector_fields, secondary_vector_fields = [], []

for vector_field_name in vector_fields:
    if field_is_primary(vector_field_name):
        primary_vector_fields.append(vector_field_name)
    elif not field_is_primary(vector_field_name):
        secondary_vector_fields.append(vector_field_name)
        
##############################

#DOWNLOAD PRIMARY (ECCO) DATA

#Iterate over primary scalar fields; download associated data
for field_name in primary_scalar_fields:
    download_data.main(field_name=field_name, initial_month=initial_month, initial_year=initial_year, final_month=final_month, final_year=final_year, time_ave_type=time_ave_type, datdir_primary=datdir_primary)
    
#Iterate over primary vector fields; download associated data
for field_name in primary_vector_fields:
    download_data.main(field_name=field_name, initial_month=initial_month, initial_year=initial_year, final_month=final_month, final_year=final_year, time_ave_type=time_ave_type, datdir_primary=datdir_primary)

##############################

#COMPUTE SECONDARY DATA

for field_name in secondary_scalar_fields: #Iterate over scalar fields; compute data if nonexistent
    create_comp_data_file(field_name, initial_month, initial_year, final_month, final_year, datdir_primary, datdir_secondary, rho_ref, nu_E, time_ave_type)

for field_name in secondary_vector_fields: #Iterate over vector fields; compute data if nonexistent
    create_comp_data_file(field_name, initial_month, initial_year, final_month, final_year, datdir_primary, datdir_secondary, rho_ref, nu_E, time_ave_type)

##############################

#ASSEMBLE A LIST OF STRINGS REPRESENTING TIMES TO ITERATE OVER

if time_ave_type == 'monthly': #will update to include other options

    date_strings = []
    month, year = int(initial_month), int(initial_year)
    
    #Append every eligible date string to list 'date_strings'
    
    while year < int(final_year):
        
        while month <= 12:
            date_string = str(year) + '-' + get_monthstr(month)
            date_strings.append(date_string)
            month += 1
            
        year += 1
        month = 1
        
    if year == int(final_year):

        while month <= int(final_month):
            date_string = final_year + '-' + get_monthstr(month)
            date_strings.append(date_string)
            month += 1
            
##############################
            
#LOAD AND VISUALIZE DATA
        
for date_string in date_strings:

    for scalar_field_name in scalar_fields:

        #Load scalar data

        if field_is_primary(scalar_field_name): #Call function that loads ECCO data
            scalar_ds = load_ECCO_data_file(scalar_field_name, date_string, datdir_primary, time_ave_type)

        elif not field_is_primary(scalar_field_name): #Call function that loads computed data
            scalar_ds = load_comp_data_file(scalar_field_name, date_string, datdir_secondary, time_ave_type)

        #Plot scalar field on its own 

        #for vector_field_name in vector_fields:

            #Plot this vector with the scalar

        #Save plots
    
##############################

#REMOVE SAVED PRIMARY DATA, IF INDICATED

#if clear_data_files: