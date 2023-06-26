"""
Saves annual averages of ECCO (or other) field, from monthly averages.

Rosalie Cormier, 2023
"""

##############################

#IMPORTS AND SETTINGS

import sys
import os
import argparse
import xarray as xr

from os.path import expanduser, join

from ecco_general import load_dataset, comp_temp_mean
from ecco_field_variables import get_field_vars

##############################

#PARSE COMMAND-LINE INPUT AND SET GLOBAL VARIABLES

parser = argparse.ArgumentParser(description="Average field over a year", \
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--field", type=str, help="Name of field to average", default="PHIHYDcR")
parser.add_argument("--datdir", type=str, help="Directory (rel. to home or here) with monthly data", \
                    default="Downloads")
parser.add_argument("--usecompdata", dest='usecompdata', help="Whether to use computed monthly data", \
                    default=False, action='store_true')
parser.add_argument("--outdir", type=str, help="Directory (rel. to here) to save output", \
                    default="yearly_averages")

parser.add_argument("years", type=int, help="Years to average over (separately)", nargs="+")

args = parser.parse_args()
config = vars(args)

field = config['field']
years = config['years']

months = ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"]

usecompdata = config['usecompdata']
    
if not usecompdata:
    
    user_home_dir = expanduser('~')
    sys.path.append(join(user_home_dir, 'ECCOv4-py'))
    datdir = join(user_home_dir, config['datdir'], 'ECCO_V4r4_PODAAC')
    
elif usecompdata:
    datdir = join(".", config['datdir'])

outdir = join(".", config['outdir'])

if not os.path.exists(outdir):
    os.makedirs(outdir)

##############################

monthly_shortname, monthly_nc_str = get_field_vars(field)
    
download_dir = join(datdir, monthly_shortname) #Get file list
    
for year in years: #Iterate over years
        
    year = str(year)
        
    monthly_fields = []
        
    for month in months: #Iterate over months
        
        if not usecompdata:
            
            curr_file = join(download_dir, monthly_nc_str+year+"-"+month+"_ECCO_V4r4_native_llc0090.nc")
            ds_month = load_dataset(curr_file) #Load monthly file into workspace
            
        elif usecompdata:
            
            curr_file = join(download_dir, monthly_nc_str+year+"-"+month+".nc")
            ds_month = xr.open_mfdataset(curr_file, engine="scipy")
            
        monthly_fields.append(ds_month.squeeze())
            
    yearly_avg_field = comp_temp_mean(monthly_fields)
    yearly_avg_field.to_netcdf(path=join(outdir, "avg_"+field+"_"+str(year)+".nc"), engine="scipy")
    yearly_avg_field.close()