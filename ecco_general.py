"""
Rosalie Cormier, 2023
"""

import xarray as xr

from os.path import join
from ecco_download import ecco_podaac_download

def load_grid(datdir):
    
    """
    Loads ECCO grid.
    
    datdir = directory where data is stored
    """
    
    grid_params_shortname = "ECCO_L4_GEOMETRY_LLC0090GRID_V4R4"

    #Download ECCO grid parameters (date is arbitrary)
    ecco_podaac_download(ShortName=grid_params_shortname, StartDate="2000-01-01", \
                     EndDate="2000-01-02", download_root_dir=datdir, n_workers=6, \
                     force_redownload=False)

    grid_params_file = "GRID_GEOMETRY_ECCO_V4r4_native_llc0090.nc"
    grid_params_file_path = join(datdir, grid_params_shortname, grid_params_file)

    #Load grid parameters
    ds_grid = xr.open_dataset(grid_params_file_path)
    
    return ds_grid

def get_month_end(monthstr, yearstr):
    
    """
    Return string representing last day of specified month.
    """
    
    month_end_dict = {"01": "31", "03": "31", "04": "30",
                 "05": "31", "06": "30", "07": "31", "08": "31",
                 "09": "30", "10": "31", "11": "30", "12": "31"}
    
    if monthstr != "02":
        endmonth = month_end_dict[monthstr]
        
    elif monthstr == "02":
        if int(yearstr) % 4 == 0:
            endmonth = "29"
        else: 
            endmonth = "28"
    
    return endmonth