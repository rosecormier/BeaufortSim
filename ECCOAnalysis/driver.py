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

from os.path import join

import read_input
import download_data
import comp_secondary
import remove_data

from functions_field_variables import get_field_lists
from functions_remove_data import get_date_strings
from functions_visualization import ArcCir_pcolormesh

##############################

#READ (AND LOG) INPUT PARAMETERS; SAVE TO VARIABLES

logfile_path, param_data = read_input.main()

#Directories
datdir_primary = param_data.data_folder_primary
datdir_secondary = param_data.data_folder_secondary
visdir = param_data.visualization_folder

#Ocean properties
rho_ref, nu_E = param_data.rho_ref, param_data.nu_E

#Temporal inputs
time_ave_type = param_data.time_ave_type
initial_month, initial_year = param_data.initial_month, param_data.initial_year
final_month, final_year = param_data.final_month, param_data.final_year
time_kwargs = param_data.time_kwargs

#Spatial inputs

plot_plane_type = param_data.plot_plane_type

if plot_plane_type == "depth_index_const":
    depth = param_data.depth_index_range[0]
    lat_res, lon_res = param_data.lat_res, param_data.lon_res
    latmin, latmax = param_data.lat_range[0], param_data.lat_range[1]
    lonmin, lonmax = param_data.lon_range[0], param_data.lon_range[1]
    plane_string = "k{}".format(depth) #To be used in filenames
    
elif plot_plane_type == "latitude_const":
    lat = param_data.lat_range[0]
    lon_res = param_data.lon_res
    depthmin = param_data.depth_index_range[0] 
    depthmax = param_data.depth_index_range[1]
    lonmin, lonmax = param_data.lon_range[0], param_data.lon_range[1]
    plane_string = "lat{}".format(lat) #To be used in filenames
    
elif plot_plane_type == "longitude_const":
    lon = param_data.lon_range[0]
    lat_res = param_data.lat_res
    depthmin = param_data.depth_index_range[0]
    depthmax = param_data.depth_index_range[1]
    latmin, latmax = param_data.lat_range[0], param_data.lat_range[1]
    plane_string = "lon{}".format(lon) #To be used in filenames
    
#Field inputs
scalar_fields = param_data.scalar_fields 
vector_fields = param_data.vector_fields

#Flag to clear primary datafiles after use
clear_data_files = param_data.clear_data_files 

##############################

#DISTINGUISH PRIMARY (ECCO) FIELDS FROM SECONDARY (COMPUTED) FIELDS

field_lists = get_field_lists(scalar_fields, vector_fields)

primary_scalar_fields = field_lists[0]
secondary_scalar_fields = field_lists[1]
primary_vector_fields = field_lists[2]
secondary_vector_fields = field_lists[3]
        
##############################

#DOWNLOAD PRIMARY (ECCO) DATA

ECCO_file_date_strings = []

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

#COMPUTE SECONDARY DATA

date_strings = get_date_strings(initial_month, initial_year, final_month, 
                                final_year, time_ave_type, time_kwargs)

for field_name in secondary_scalar_fields: #Iterate over scalar fields
    for date_string in date_strings: 
        #Compute data if nonexistent
        comp_secondary.main(datdir_primary=datdir_primary, 
                            datdir_secondary=datdir_secondary, 
                            date_string=date_string, 
                            time_ave_type=time_ave_type, 
                            time_kwargs=time_kwargs, field_name=field_name, 
                            rho_ref=rho_ref, nu_E=nu_E, depth_index=int(depth))
        
for field_name in secondary_vector_fields: #Iterate over vector fields
    for date_string in date_strings:
        #Compute data if nonexistent
        comp_secondary.main(datdir_primary=datdir_primary, 
                            datdir_secondary=datdir_secondary, 
                            date_string=date_string, 
                            time_ave_type=time_ave_type, 
                            time_kwargs=time_kwargs, field_name=field_name, 
                            rho_ref=rho_ref, nu_E=nu_E, depth_index=int(depth))

##############################
        
#VISUALIZE DATA

spatial_bounds = [depth, latmin, latmax, lonmin, lonmax] 
resolutions = [lat_res, lon_res]
#Will update these lines when I modify to allow other plane types
        
for date_string in date_strings:

    for scalar_field_name in scalar_fields:
        
        #Set up output directory and file to save scalar plot to
        outfile = join(visdir, plot_plane_type, time_ave_type, 
                       '{}_{}_{}.png'.format(scalar_field_name, plane_string, 
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
                           '{}_{}_{}_{}.png'.format(scalar_field_name, 
                                                    vector_field_name, 
                                                    plane_string, date_string))
            
            #Plot this vector with the scalar
            ArcCir_pcolormesh(scalar_field_name, date_string, datdir_primary, 
                              datdir_secondary, time_ave_type, plot_plane_type, 
                              spatial_bounds, resolutions, outfile, 
                              vector_field_name=vector_field_name)
    
##############################

#REMOVE PRIMARY DATAFILES, IF INDICATED

if clear_data_files:

    if time_ave_type == 'seasonal': 
        #In this case, data will have been downloaded for every month in each 
        #season, so update the set of date strings to loop through
        date_strings = ECCO_file_date_strings 

    remove_data.main(clear_data_files=clear_data_files, date_strings=date_strings, 
                     primary_scalar_fields=primary_scalar_fields, 
                     primary_vector_fields=primary_vector_fields, 
                     datdir_primary=datdir_primary, 
                     secondary_scalar_fields=secondary_scalar_fields, 
                     secondary_vector_fields=secondary_vector_fields)