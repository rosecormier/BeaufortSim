"""
Rosalie Cormier, 2023
"""

##############################

#IMPORTS AND SETTINGS

import sys
import os
import argparse
import ecco_v4_py as ecco
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

while i < mos:
    
    month = month_dict[i % 12]
    endmonth = month_end_dict[month]
    
    #Download the monthly-averaged velocity file
    ecco_podaac_download(ShortName=vel_monthly_shortname, StartDate="2000-"+month+"-02", 
                         EndDate="2000-"+month+"-"+endmonth, download_root_dir=datdir, n_workers=6, 
                         force_redownload=False)
    
    #Download the monthly-averaged density-/pressure-anomaly file
    ecco_podaac_download(ShortName=denspress_monthly_shortname, StartDate="2000-"+month+"-02", 
                         EndDate="2000-"+month+"-"+endmonth, download_root_dir=datdir, n_workers=6, 
                         force_redownload=False)
    
    i += 1 #Next month
    
#Download ECCO grid parameters
ecco_podaac_download(ShortName=grid_params_shortname, StartDate="2000-01-01", EndDate="2000-01-01", 
                     download_root_dir=datdir, n_workers=6, force_redownload=False)