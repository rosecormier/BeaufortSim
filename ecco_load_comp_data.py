"""
Loads a specified NetCDF file containing a DataSet computed from ECCO data.

Rosalie Cormier, 2023
"""

import os
import xarray as xr

from ecco_general import load_dataset

import compute_monthly_avgs

##############################

def load_comp_file(monthly_file, latmin, latmax, lonmin, lonmax, year, datdir, compdatdir):
    
    """
    Checks that a computed file exists, and creates it if it doesn't, then loads DataSet.
    """
    
    if os.path.exists(monthly_file): #Look for the file
        ds_month = xr.open_mfdataset(monthly_file, engine="scipy") #Load DataSet

    else: #If it doesn't exist, compute it (for the year)
                                
        compute_monthly_avgs.main(latmin=latmin, latmax=latmax, lonmin=lonmin, lonmax=lonmax, startyr=year, years=1, datdir=datdir, outdir=compdatdir)
        
        ds_month = load_dataset(monthly_file) #Load DataSet
        
    return ds_month