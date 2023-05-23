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

from os.path import expanduser, join
from ecco_download import ecco_podaac_download

from ecco_visualization import *

##############################

#PARSE COMMAND-LINE INPUT AND SET GLOBAL VARIABLES

parser = argparse.ArgumentParser(description="Plot pressure and velocity fields in Beaufort Gyre",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--lats", type=float, help="Bounding latitudes", nargs=2, default=[70.0, 85.0])
parser.add_argument("--lons", type=float, help="Bounding longitudes", nargs=2, default=[-180.0, -90.0])
parser.add_argument("--month", type=str, help="Start month", default="01")
parser.add_argument("--months", type=int, help="Total number of months", default=12)
parser.add_argument("--kvals", type=int, help="Bounding k-values", nargs=2, default=[0, 4])
parser.add_argument("--res", type=float, help="Lat/lon resolution in degrees", nargs=1, default=1.0)
parser.add_argument("--datdir", type=str, help="Directory (relative to home) to store ECCO data", default="Downloads")
parser.add_argument("--outdir", type=str, help="Output directory (relative to here)", default="visualization")

parser.add_argument("start", type=int, help="Start year")

args = parser.parse_args()

config = vars(args)

latmin, latmax = config['lats'][0], config['lats'][1]
lonmin, lonmax = config['lons'][0], config['lons'][1]
startmo, startyr, mos = config['month'], config['start'], config['months']
kmin, kmax = config['kvals'][0], config['kvals'][1]
resolution = config['res']

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
monthstrs, yearstrs = [], []

#Iterate over all specified months
while i < mos:

    monthstr = month_dict[i % 12]
    endmonth = month_end_dict[monthstr]
    yearstr = str(year)
    
    monthstrs.append(monthstr)
    yearstrs.append(yearstr)
    
    #Download the monthly-averaged velocity file
    ecco_podaac_download(ShortName=vel_monthly_shortname, StartDate=yearstr+"-"+monthstr+"-02", 
                         EndDate=yearstr+"-"+monthstr+"-"+endmonth, download_root_dir=datdir, n_workers=6, 
                         force_redownload=False)
     
    #Download the monthly-averaged density-/pressure-anomaly file
    ecco_podaac_download(ShortName=denspress_monthly_shortname, StartDate=yearstr+"-"+monthstr+"-02", 
                         EndDate=yearstr+"-"+monthstr+"-"+endmonth, download_root_dir=datdir, n_workers=6, 
                         force_redownload=False)
    
    if i == "12":
        year += 1 #Go to next year
        
    i += 1 #Go to next month
    
#Download ECCO grid parameters (date is arbitrary)
ecco_podaac_download(ShortName=grid_params_shortname, StartDate="2000-01-01", EndDate="2000-01-02", 
                     download_root_dir=datdir, n_workers=6, force_redownload=False)

##############################

#GET FILE LISTS AND LOAD GRID

vel_dir = join(datdir, vel_monthly_shortname)
denspress_dir = join(datdir, denspress_monthly_shortname)

grid_params_file = "GRID_GEOMETRY_ECCO_V4r4_native_llc0090.nc"
grid_params_file_path = join(datdir, grid_params_shortname, grid_params_file)

#Load grid parameters
ds_grid = xr.open_dataset(grid_params_file_path)

##############################

#CREATE MONTHLY PLOTS OF VELOCITY AND PRESSURE ANOMALY

vir_nanmasked = plt.get_cmap('viridis_r').copy()
vir_nanmasked.set_bad('black')

ds_vels, ds_pressures = [], []

for m in range(mos):
    
    monthstr = monthstrs[m]
    yearstr = yearstrs[m]
    
    curr_vel_file = join(vel_dir, "OCEAN_VELOCITY_mon_mean_"+yearstr+
                         "-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
    curr_denspress_file = join(denspress_dir, "OCEAN_DENS_STRAT_PRESS_mon_mean_"+
                               yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")

    #Load monthly velocity file into workspace
    ds_vel_mo = xr.open_mfdataset(curr_vel_file, parallel=True, data_vars='minimal', coords='minimal', 
                              compat='override')
    
    ds_vels.append(ds_vel_mo) 
    
    #Interpolate velocities to centres of grid cells
    ds_vel_mo.UVEL.data, ds_vel_mo.VVEL.data = ds_vel_mo.UVEL.values, ds_vel_mo.VVEL.values
    
    #Load monthly density-/pressure-anomaly file into workspace
    ds_denspress_mo = xr.open_mfdataset(curr_denspress_file, parallel=True, data_vars='minimal', coords='minimal', 
                                    compat='override')
    
    ds_pressures.append(ds_denspress_mo)

    #Plot velocity and pressure fields
    
    press_vel_title = 'Pressure anomaly and water velocity in Arctic Circle \n at {} ({}-{})'.format('replace', \
                                                                                                yearstr, \
                                                                                                monthstr)
    
    ArcCir_contourf_quiver(ds_grid, 1, [ds_denspress_mo], [ds_vel_mo], 'PHIHYDcR', 'UVEL', 'VVEL', resolution, vir_nanmasked, [93, 97], yearstr+"-"+monthstr, outfile=join(outdir, 'u_p_anom_{}-{}.pdf'.format(monthstr, yearstr)), latmin=latmin, latmax=latmax, lonmin=lonmin, lonmax=lonmax)

#Plot all months
ArcCir_contourf_quiver_grid(ds_grid, 1, ds_pressures, ds_vels, 'PHIHYDcR', [93, 97], 'UVEL', 'VVEL', resolution, 
                           vir_nanmasked, monthstrs, yearstrs, outfile=join(outdir, 'u_p_anom_all{}.png'.format(yearstr)),
                           latmin=latmin, latmax=latmax, lonmin=lonmin, lonmax=lonmax)

#Plot annual averages
ArcCir_contourf_quiver(ds_grid, 1, ds_pressures, ds_vels, 'PHIHYDcR', 'UVEL', 'VVEL', resolution, vir_nanmasked, [93, 97], yearstrs[0]+" average", outfile=join(outdir, 'u_p_anom_avg{}.pdf'.format(yearstr)), latmin=latmin, latmax=latmax, lonmin=lonmin, lonmax=lonmax)
"""   
#Compute and plot residuals of annual averages

ds_press_residuals = []
ds_vel_residuals = []

for pressure in ds_pressures:
    residual = pressure.copy() * 0
    concat_pressure = xr.concat((pressure, -1*ds_press_mean), dim='time')
    residual_press = concat_pressure['PHIHYDcR'].sum(dim=['time'])
    residual['PHIHYDcR'] = residual_press
    ds_press_residuals.append(residual)
    
for vel in ds_vels:
    residual = vel.copy() * 0
    concat_vel = xr.concat((vel, -1*ds_vel_mean), dim='time')
    residual_u = concat_vel['UVEL'].sum(dim=['time'])
    residual_v = concat_vel['VVEL'].sum(dim=['time'])
    residual['UVEL'] = residual_u
    residual['VVEL'] = residual_v
    ds_vel_residuals.append(residual)

#Title for residual plot
title = ''

#Plot residuals (pressure and velocity) for all months
ArcCir_contourf_quiver_grid(ds_grid, 1, ds_press_residuals, ds_vel_residuals, 'PHIHYDcR', [-2, 2], 'UVEL', 'VVEL', resolution, 
                           vir_nanmasked, monthstrs, yearstrs, outfile=join(outdir, 'u_p_resids_all{}.png'.format(yearstr)),
                           latmin=latmin, latmax=latmax, lonmin=lonmin, lonmax=lonmax)
                           
"""