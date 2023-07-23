"""
Load monthly data from a previously downloaded ECCO dataset.

Rosalie Cormier, 2023
"""

import os

from os.path import join

from ecco_general import load_dataset

import download_new_data

##############################

def load_ECCO_dataset(variable_dir, variable_monthly_nc_str, yearstr, monthstr, year, datdir, scalar_attr, xvec_attr):
    
    curr_file = join(variable_dir, variable_monthly_nc_str+yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")

    if not os.path.exists(curr_file): #If the file doesn't exist, download it
                
        if scalar_attr is not None:
            download_new_data.main(startmo="01", startyr=year, months=12, scalars=[scalar_attr], xvectors=None, datdir=datdir)
        elif xvec_attr is not None:
            download_new_data.main(startmo="01", startyr=year, months=12, scalars=None, xvectors=[xvec_attr], datdir=datdir)
        
    ds_month = load_dataset(curr_file)
    
    return ds_month

##############################

def main(**kwargs):
    
    #Define variables
    
    variable_dir = kwargs.get('variable_dir')
    variable_monthly_nc_str = kwargs.get('variable_monthly_nc_str')
    yearstr = kwargs.get('yearstr')
    monthstr = kwargs.get('monthstr')
    year = kwargs.get('year')
    datdir = kwargs.get('datdir')
    scalar_attr = kwargs.get('scalar_attr')
    xvec_attr = kwargs.get('xvec_attr')
    
    load_ECCO_dataset(variable_dir, variable_monthly_nc_str, yearstr, monthstr, year, datdir, scalar_attr=scalar_attr, xvec_attr=xvec_attr)
    
##############################

if __name__ == "__main__":
    main()