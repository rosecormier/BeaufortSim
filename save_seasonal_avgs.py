"""
Saves seasonal averages of ECCO fields, from monthly averages.

Rosalie Cormier, 2023
"""

##############################

#IMPORTS AND SETTINGS

import sys
import os
import argparse

from os.path import expanduser, join

from ecco_general import load_dataset, comp_temp_mean, get_season_months_and_years
from ecco_field_variables_new import get_field_vars

##############################

#PARSE COMMAND-LINE INPUT AND SET GLOBAL VARIABLES

parser = argparse.ArgumentParser(description="Average ECCO fields over a season", \
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--ECCOfields", type=str, help="Names of ECCO fields to average", nargs="+")
parser.add_argument("--months", type=str, help="Months marking start/end of season", nargs=2, \
                    default=["01", "01"])
parser.add_argument("--datdir", type=str, help="Directory (rel. to home) to store ECCO data", \
                    default="Downloads")
parser.add_argument("--outdir", type=str, help="Directory (rel. to here) to save output", \
                    default="seasonal_averages")

parser.add_argument("years", type=int, help="Years to average over (separately)", nargs="+")

args = parser.parse_args()
config = vars(args)

ECCO_fields = config['ECCOfields']
years = config['years']
start_month, end_month = config['months']

user_home_dir = expanduser('~')
sys.path.append(join(user_home_dir, 'ECCOv4-py'))
datdir = join(user_home_dir, config['datdir'], 'ECCO_V4r4_PODAAC')

outdir = join(".", config['outdir'])

if not os.path.exists(outdir):
    os.makedirs(outdir)

##############################

for field in ECCO_fields: #Iterate over fields
    
    monthly_shortname, monthly_nc_str = get_field_vars(field)
    
    download_dir = join(datdir, monthly_shortname) #Get file list
    
    for i in range(len(years)): #Iterate over seasons
        
        year_start = years[i]
        season_months, season_years = get_season_months_and_years(start_month, end_month)
        
        monthly_fields = []

        for j in range(len(season_months)):
            
            month, year_actual = season_months[j], year_start + season_years[j]
            year = str(year_actual)
            
            curr_file = join(download_dir, monthly_nc_str+year+"-"+month+"_ECCO_V4r4_native_llc0090.nc")
            
            ds_month = load_dataset(curr_file) #Load monthly file into workspace
            
            monthly_fields.append(ds_month.squeeze())
        
        seasonal_avg_field = comp_temp_mean(monthly_fields)
        seasonal_avg_field.to_netcdf(path=join(outdir, "avg_"+field+"_"+start_month+str(year_start)+"-"+end_month+year+".nc"), engine="scipy")
        seasonal_avg_field.close()