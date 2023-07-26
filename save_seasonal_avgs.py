"""
Saves seasonal averages of ECCO (or other) field, from monthly averages.

Rosalie Cormier, 2023
"""

##############################

#IMPORTS AND SETTINGS

import sys
import os
import argparse
import xarray as xr

from os.path import expanduser, join

from functions_ecco_general import load_dataset, comp_temp_mean, get_season_months_and_years
from functions_field_variables import get_field_vars

#To be called from this script

import compute_monthly_avgs
import download_new_data

##############################

def main(**kwargs):
    
    if not kwargs:

        parser = argparse.ArgumentParser(description="Average ECCO fields over a season", \
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        parser.add_argument("--field", type=str, help="Name of field to average", default="PHIHYDcR")
        parser.add_argument("--months", type=str, help="Months marking start/end of season", nargs=2, \
                            default=["01", "01"])
        parser.add_argument("--datdir", type=str, help="Directory (rel. to home or here) with monthly data", \
                            default="Downloads")
        parser.add_argument("--usecompdata", dest='usecompdata', help="Whether to use computed monthly data", \
                            default=False, action='store_true')
        parser.add_argument("--outdir", type=str, help="Directory (rel. to here) to save output", \
                            default="seasonal_averages")

        parser.add_argument("years", type=int, help="Years to average over (separately)", nargs="+")

        args = parser.parse_args()
        config = vars(args)

        field = config['field']
        years = config['years']
        start_month, end_month = config['months']

        usecompdata = config['usecompdata']
        
        datdirshort = config['datdir']
        outdir = config['outdir']
        
    elif kwargs:
        
        field = kwargs.get('field')
        years = kwargs.get('years')
        start_month, end_month = kwargs.get('start_month'), kwargs.get('end_month')
        usecompdata = kwargs.get('usecompdata')
        datdirshort = kwargs.get('datdir')
        outdir = kwargs.get('outdir')

    if not usecompdata:

        homedir = expanduser('~')
        sys.path.append(join(homedir, 'ECCOv4-py'))
        datdir = join(homedir, datdirshort, 'ECCO_V4r4_PODAAC')

    elif usecompdata:
        datdir = join(".", datdirshort)

    outdir = join(".", outdir)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    ##############################

    monthly_shortname, monthly_nc_str = get_field_vars(field)

    download_dir = join(datdir, monthly_shortname) #Get file list

    for i in range(len(years)): #Iterate over seasons

        year_start = years[i]
        season_months, season_years = get_season_months_and_years(start_month, end_month)

        monthly_fields = []

        for j in range(len(season_months)):

            month, year_actual = season_months[j], year_start + season_years[j]
            year = str(year_actual)

            if not usecompdata:

                curr_file = join(download_dir, monthly_nc_str+year+"-"+month+"_ECCO_V4r4_native_llc0090.nc")
                
                if not os.path.exists(curr_file): #Download if it doesn't exist
                    download_new_data.main(startmo=start_month, startyr=year_start, months=len(season_months), scalars=[field], xvectors=None, datdir=datdirshort)
                
                ds_month = load_dataset(curr_file) #Load monthly file into workspace

            elif usecompdata:

                curr_file = join(download_dir, monthly_nc_str+year+"-"+month+".nc")
                
                if not os.path.exists(curr_file): #Compute if it doesn't exist
                    compute_monthly_avgs.main(latmin=70.0, latmax=85.0, lonmin=-180.0, lonmax=-90.0, startyr=year_start, years=2, datdir='Downloads', outdir='computed_monthly')
                
                ds_month = xr.open_mfdataset(curr_file, engine="scipy")

            monthly_fields.append(ds_month.squeeze())

        seasonal_avg_field = comp_temp_mean(monthly_fields)
        seasonal_avg_field.to_netcdf(path=join(outdir, "avg_"+field+"_"+start_month+str(year_start)+"-"+end_month+year+".nc"), engine="scipy")
        seasonal_avg_field.close()
        
##############################

if __name__ == "__main__":
    main()