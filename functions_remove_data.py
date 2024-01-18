"""
Function to remove primary data.

Rosalie Cormier, 2024
"""

import os
import glob

from os.path import join

from functions_field_variables import get_field_variable, \
get_monthly_shortname, get_monthly_nc_string, get_seasonal_shortname, \
get_seasonal_nc_string

##############################

def remove_primary_files(field_name, datdir_primary, date_strings):

    field_shortname = get_monthly_shortname(get_field_variable(field_name))
    field_nc_string = get_monthly_nc_string(get_field_variable(field_name))

    #To avoid errors, only remove files after confirming directory exists
    if os.path.exists(join(datdir_primary, field_shortname)):
        for date_string in date_strings: #Iterate over times
            pattern = join(datdir_primary, field_shortname, 
                           field_nc_string+date_string+r"*")
            for item in glob.iglob(pattern, recursive=True): 
                os.remove(item) #Delete the files
        try:
            #Delete the directory, if empty
            os.rmdir(join(datdir_primary, field_shortname)) 
        except:
            return join(datdir_primary, field_shortname)