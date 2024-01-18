"""
This script uses the data provided in the input file (input.txt) to do the 
following:
    -Save a txt file listing the input parameters to be used;
    -Download data corresponding to specified fields;
    -Compute and save secondary fields (skips if unnecessary or done);
    -Plot data for all fields specified and save plots; and
    -Remove primary data files.

R. Cormier, F. Poulin, 2024
"""

import os
import sys

from os.path import expanduser, join

import read_input
import download_data
import comp_secondary

from functions_ecco_general import get_monthstr
from functions_field_variables import field_is_primary
from functions_remove_data import remove_primary_files
from functions_visualization import ArcCir_pcolormesh

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

experiment_number = 0 
#Hardcoded for now; eventually will need to automatically produce a unique 
#number for each run

for directory in [visdir, logdir]:
    if not os.path.exists(directory):
        os.makedirs(directory)

#Create the log file
logfile = join(logdir, "logfile_{}.txt".format(str(experiment_number))) 

f = open(logfile, "w") 
#When ready to use the script in full, switch "w" to "x"

#Ocean properties

rho_ref = param_data.rho_ref
nu_E = param_data.nu_E
f.write("rho_ref (kg/m^3) = " + rho_ref + "\n" + "nu_E (m^2/s) = " + nu_E + 
        "\n\n")

#Temporal inputs

time_ave_type = param_data.time_ave_type
initial_year = param_data.initial_month[0]
initial_month = param_data.initial_month[1]
final_year, final_month = param_data.final_month[0], param_data.final_month[1]
f.write("time_ave_type = " + time_ave_type + "\n")
f.write("initial_month = " + initial_month + "\n" + "initial_year = " + 
        initial_year + "\n")
f.write("final_month = " + final_month + "\n" + "final_year = " + final_year 
        + "\n\n")

if time_ave_type == "daily":
    initial_day, final_day = param_data.initial_day, param_data.final_day
    f.write("initial_day = " + initial_day + "\n" + "final_day = " + 
            final_day + "\n\n")
    time_kwargs = [initial_day, final_day]
elif time_ave_type == "seasonal":
    season_start, season_end = param_data.season_start, param_data.season_end
    f.write("season_start = " + season_start + "\n" + "season_end = " + 
            season_end + "\n\n")
    time_kwargs = [season_start, season_end]
else:
    time_kwargs = None

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
    
    plane_string = 'k'+depth #To be used in filenames
    
elif plot_plane_type == "latitude_const":
    
    lat = param_data.lat_range[0]
    lon_res = param_data.lon_res
    depthmin, depthmax = param_data.depth_range[0], param_data.depth_range[1]
    lonmin, lonmax = param_data.lon_range[0], param_data.lon_range[1]
    f.write("lat = " + lat + "\n")
    f.write("lon_res = " + lon_res + "\n")
    f.write("depthmin = " + depthmin + "\n" + "depthmax = " + depthmax + "\n")
    f.write("lonmin = " + lonmin + "\n" + "lonmax = " + lonmax + "\n")
    
    plane_string = 'lat'+lat #To be used in filenames
    
elif plot_plane_type == "longitude_const":
    
    lon = param_data.lon_range[0]
    lat_res = param_data.lat_res
    depthmin, depthmax = param_data.depth_range[0], param_data.depth_range[1]
    latmin, latmax = param_data.lat_range[0], param_data.lat_range[1]
    f.write("lon = " + lon + "\n")
    f.write("lat_res = " + lat_res + "\n")
    f.write("depthmin = " + depthmin + "\n" + "depthmax = " + depthmax + "\n")
    f.write("latmin = " + latmin + "\n" + "latmax = " + latmax + "\n")
    
    plane_string = 'lon'+lon #To be used in filenames
    
# interp_type = #TBA

#Field inputs

scalar_fields = param_data.scalar_fields
f.write("scalar_fields = ")
for i in range(len(scalar_fields)-1):
    f.write(scalar_fields[i] + ", ")
f.write(scalar_fields[-1] + "\n")
        
vector_fields = param_data.vector_fields
f.write("vector_fields = ")
if vector_fields == ['']:
    vector_fields = []
else:
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

if len(vector_fields) != 0:
    for vector_field_name in vector_fields:
        if field_is_primary(vector_field_name):
            primary_vector_fields.append(vector_field_name)
        elif not field_is_primary(vector_field_name):
            secondary_vector_fields.append(vector_field_name)
        
##############################

#DOWNLOAD PRIMARY (ECCO) DATA

#Iterate over primary scalar fields; download associated data
for field_name in primary_scalar_fields:
    ECCO_file_date_strings = download_data.main(field_name=field_name, 
                                                initial_month=initial_month, 
                                                initial_year=initial_year, 
                                                final_month=final_month, 
                                                final_year=final_year, 
                                                time_ave_type=time_ave_type, 
                                                datdir_primary=datdir_primary, 
                                                time_kwargs=time_kwargs)
    
#Iterate over primary vector fields; download associated data
for field_name in primary_vector_fields:
    ECCO_file_date_strings = download_data.main(field_name=field_name, 
                                                initial_month=initial_month, 
                                                initial_year=initial_year, 
                                                final_month=final_month, 
                                                final_year=final_year, 
                                                time_ave_type=time_ave_type, 
                                                datdir_primary=datdir_primary, 
                                                time_kwargs=time_kwargs)

##############################

#ASSEMBLE A LIST OF STRINGS REPRESENTING TIMES TO ITERATE OVER

date_strings = []
month, year = int(initial_month), int(initial_year)

#Append every eligible date string to list 'date_strings'

if time_ave_type == 'monthly':
    
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
            
elif time_ave_type == 'seasonal':
        
    if int(season_start) < int(season_end):
        while year <= int(final_year):
            date_string = '{}-{}_{}'.format(season_start, season_end, str(year))
            date_strings.append(date_string) 
            year += 1
            
    elif int(season_end) < int(season_start):
        while year < int(final_year):
            date_string = '{}-{}_{}-{}'.format(season_start, season_end, 
                                               str(year), str(year+1))
            date_strings.append(date_string)
            year += 1

##############################

#COMPUTE SECONDARY DATA

for field_name in secondary_scalar_fields: #Iterate over scalar fields
    for date_string in date_strings: #Iterate over times 
        #Compute data if nonexistent
        comp_secondary.main(datdir_primary=datdir_primary, 
                            datdir_secondary=datdir_secondary, 
                            date_string=date_string, 
                            time_ave_type=time_ave_type, 
                            time_kwargs=time_kwargs, field_name=field_name, 
                            rho_ref=rho_ref, nu_E=nu_E)
        
for field_name in secondary_vector_fields: #Iterate over vector fields
    for date_string in date_strings: #Iterate over times
        #Compute data if nonexistent
        comp_secondary.main(datdir_primary=datdir_primary, 
                            datdir_secondary=datdir_secondary, 
                            date_string=date_string, 
                            time_ave_type=time_ave_type, 
                            time_kwargs=time_kwargs, field_name=field_name, 
                            rho_ref=rho_ref, nu_E=nu_E)

##############################
            
#VISUALIZE DATA

spatial_bounds = [depth, latmin, latmax, lonmin, lonmax] 
resolutions = [lat_res, lon_res]
#Will update these lines when I modify to allow other plane types
        
for date_string in date_strings: #Iterate over times

    for scalar_field_name in scalar_fields:
        
        #Set up output directory and file to save scalar plot to
        
        outfile = join(visdir, plot_plane_type, time_ave_type, 
                       '{}_{}_{}.pdf'.format(scalar_field_name, plane_string, 
                                             date_string))
        if not os.path.exists(join(visdir, plot_plane_type, time_ave_type)):
            os.makedirs(join(visdir, plot_plane_type, time_ave_type))
        
        #Plot scalar field on its own
        ArcCir_pcolormesh(scalar_field_name, date_string, datdir_primary, 
                          datdir_secondary, time_ave_type, plot_plane_type, 
                          spatial_bounds, resolutions, outfile)

        for vector_field_name in vector_fields:

            #Set up file to save scalar/vector plot to
            outfile = join(visdir, plot_plane_type, time_ave_type, 
                           '{}_{}_{}_{}.pdf'.format(scalar_field_name, 
                                                    vector_field_name, 
                                                    plane_string, date_string))
            
            #Plot this vector with the scalar
            ArcCir_pcolormesh(scalar_field_name, date_string, datdir_primary, 
                              datdir_secondary, time_ave_type, plot_plane_type, 
                              spatial_bounds, resolutions, outfile, 
                              vector_field_name=vector_field_name)
    
##############################

#REMOVE SAVED PRIMARY DATA, IF INDICATED

if clear_data_files:
    
    if time_ave_type == 'seasonal': 
        #In this case, data will have been downloaded for every month in each 
        #season
        
        #Update the set of date strings to loop through
        date_strings = ECCO_file_date_strings 
        
    for scalar_field_name in primary_scalar_fields: #Delete primary scalar data
        remove_primary_files(scalar_field_name, datdir_primary, date_strings)
        
    print("Deleted scalar data.")
    
    for vector_field_name in primary_vector_fields: #Delete primary vector data
        remove_primary_files(vector_field_name, datdir_primary, date_strings)
        
    print("Deleted vector data.")
    
    #Identify any secondary fields that required primary data to be saved, and 
    #delete the primary data
    
    for scalar_field_name in secondary_scalar_fields:
        if scalar_field_name in ['vorticity', 'normal_strain', 'shear_strain', 
                                 '2D_div_vel']:
            remove_primary_files('horizontal_vel', datdir_primary, date_strings)
        
    for vector_field_name in secondary_vector_fields:
        if vector_field_name in ['geostrophic_vel', 'Ek_vel']:
            remove_primary_files('density_anom', datdir_primary, date_strings)
        if vector_field_name in ['Ek_vel']:
            remove_primary_files('wind_stress', datdir_primary, date_strings)
                    
    print("Done deleting primary data.")