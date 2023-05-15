"""
Rosalie Cormier, 2023
"""

##############################

#IMPORTS AND SETTINGS

import sys
import os
import argparse
import ecco_v4_py as ecco
import glob
import xarray as xr
import matplotlib.pyplot as plt

from os.path import expanduser, join
from ecco_download import ecco_podaac_download

from ecco_visualization import *

plt.rcParams["font.size"] = 12

##############################

#PARSE COMMAND-LINE INPUT AND SET GLOBAL VARIABLES

parser = argparse.ArgumentParser(description="Plot pressure and velocity fields in Beaufort Gyre",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--lats", type=float, help="Bounding latitudes", nargs=2, default=[65.0, 85.0])
parser.add_argument("--lons", type=float, help="Bounding longitudes", nargs=2, default=[90.0, 180.0])
parser.add_argument("--month", type=str, help="Start month", default="01")
parser.add_argument("--months", type=int, help="Total number of months", default=12)
parser.add_argument("--kvals", type=int, help="Bounding k-values", nargs=2, default=[0, 4])
parser.add_argument("--datdir", type=str, help="Directory (relative to home) to store ECCO data", default="Downloads")
parser.add_argument("--outdir", type=str, help="Output directory (relative to here)", default="visualization")

parser.add_argument("start", type=int, help="Start year")

args = parser.parse_args()

config = vars(args)

latmin, latmax = config['lats'][0], config['lats'][1]
lonmin, lonmax = config['lons'][0], config['lons'][1]
startmo, startyr, mos = config['month'], config['start'], config['months']
kmin, kmax = config['kvals'][0], config['kvals'][1]

user_home_dir = expanduser('~')
sys.path.append(join(user_home_dir, 'ECCOv4-py'))

datdir = join(user_home_dir, config['datdir'], 'ECCO_V4r4_PODAAC')

if not os.path.exists(datdir):
    os.makedirs(datdir)
    
outdir = join(".", config['outdir'])

if not os.path.exists(outdir):
    os.makedirs(outdir)
    
month_dict = {0: "01", 1: "02", 2: "03", 3: "04", 4: "05", 5: "06",
              6: "07", 7: "08", 8: "09", 9: "10", 10: "11", 11: "12"}
month_key_list = list(month_dict.keys())
month_val_list = list(month_dict.values())
month_end_dict = {"01": "31", "02": "28", "03": "31", "04": "30",
                 "05": "31", "06": "30", "07": "31", "08": "31",
                 "09": "30", "10": "31", "11": "30", "12": "31"} #Ideally fix Feb. at some point
    
##############################

#DOWNLOAD ECCO FILES FOR SPECIFIED MONTHS

vel_monthly_shortname = "ECCO_L4_OCEAN_VEL_LLC0090GRID_MONTHLY_V4R4"
denspress_monthly_shortname = "ECCO_L4_DENS_STRAT_PRESS_LLC0090GRID_MONTHLY_V4R4"
grid_params_shortname = "ECCO_L4_GEOMETRY_LLC0090GRID_V4R4"

month_index = month_val_list.index(startmo)
i = month_key_list[month_index]
year = startyr

#Iterate over all specified months
while i < mos:

    month = month_dict[i % 12]
    endmonth = month_end_dict[month]
    yearstr = str(year)
    
    #Download the monthly-averaged velocity file
    ecco_podaac_download(ShortName=vel_monthly_shortname, StartDate=yearstr+"-"+month+"-02", 
                         EndDate=yearstr+"-"+month+"-"+endmonth, download_root_dir=datdir, n_workers=6, 
                         force_redownload=False)
    
    #Download the monthly-averaged density-/pressure-anomaly file
    ecco_podaac_download(ShortName=denspress_monthly_shortname, StartDate=yearstr+"-"+month+"-02", 
                         EndDate=yearstr+"-"+month+"-"+endmonth, download_root_dir=datdir, n_workers=6, 
                         force_redownload=False)
    
    if i == "12":
        year += 1 #Go to next year
        
    i += 1 #Go to next month
    
#Save final date
endmo, endyr = month_dict[i % 12], yearstr
    
#Download ECCO grid parameters
ecco_podaac_download(ShortName=grid_params_shortname, StartDate="2000-01-01", EndDate="2000-01-01", 
                     download_root_dir=datdir, n_workers=6, force_redownload=False)

##############################

#LOAD DOWNLOADED FILES

vel_dir = join(datdir, vel_monthly_shortname)
curr_vel_files = list(glob.glob(join(vel_dir, '*nc')))

#Load velocity file into workspace
ds_vel_mo = xr.open_mfdataset(curr_vel_files, parallel=True, data_vars='minimal', coords='minimal', 
                              compat='override')

denspress_dir = join(datdir, denspress_monthly_shortname)
curr_denspress_files = list(glob.glob(join(denspress_dir, '*nc')))

#Load density-/pressure-anomaly file into workspace
ds_denspress_mo = xr.open_mfdataset(curr_denspress_files, parallel=True, data_vars='minimal', coords='minimal', 
                                    compat='override')

grid_params_file = "GRID_GEOMETRY_ECCO_V4r4_native_llc0090.nc"
grid_params_file_path = join(datdir, grid_params_shortname, grid_params_file)

#Load grid parameters
ds_grid = xr.open_dataset(grid_params_file_path)