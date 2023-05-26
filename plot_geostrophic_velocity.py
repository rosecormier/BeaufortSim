"""
Rosalie Cormier, 2023, based on code by Andrew Delman
"""

##############################

#IMPORTS AND SETTINGS
 
import xgcm
import os
import sys
import argparse
import xarray as xr
import matplotlib.pyplot as plt #
import matplotlib.colors as colors #
import ecco_v4_py as ecco #
import cartopy.crs as ccrs #
import numpy as np #
import cartopy.feature as cfeature #

from os.path import expanduser, join
from ecco_download import ecco_podaac_download

from ecco_general import load_grid, get_monthstr, get_month_end, get_starting_i, load_dataset
from ecco_field_variables import get_scalar_field_vars, get_vector_field_vars
from geostrophic_functions import comp_geos_vel, comp_delta_u_norm

plt.rcParams['font.size'] = 12
plt.rcParams['text.usetex'] = True

##############################

#PARSE COMMAND-LINE INPUT AND SET GLOBAL VARIABLES

parser = argparse.ArgumentParser(description="Plot geostrophic velocity in Beaufort Gyre", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--month", type=str, help="Start month", default="01")
parser.add_argument("--months", type=int, help="Total number of months", default=12)
parser.add_argument("--kvals", type=int, help="Bounding k-values", nargs=2, default=[20, 21])
parser.add_argument("--datdir", type=str, help="Directory (rel. to home) to store ECCO data", default="Downloads")
parser.add_argument("--outdir", type=str, help="Output directory (rel. to here)", default="visualization")

parser.add_argument("start", type=int, help="Start year")

args = parser.parse_args()
config = vars(args)

startmo, startyr, mos = config['month'], config['start'], config['months']
kmin, kmax = config['kvals'][0], config['kvals'][1]

user_home_dir = expanduser('~')
sys.path.append(join(user_home_dir, 'ECCOv4-py'))
datdir = join(user_home_dir, config['datdir'], 'ECCO_V4r4_PODAAC')

if not os.path.exists(datdir):
    os.makedirs(datdir)
    
outdir = join(".", config['outdir'], 'geostrophic')

if not os.path.exists(outdir):
    os.makedirs(outdir)
    
#Get variables associated with velocity and density/pressure

vel_monthly_shortname, vel_monthly_nc_str, vel_variable = get_vector_field_vars('UVEL', 'VVEL')
denspress_monthly_shortname, denspress_monthly_nc_str, denspress_variable = get_scalar_field_vars('PHIHYDcR')

rho_ref = 1029.0 #Reference density (kg/m^3)

##############################

#LOAD GRID AND DOWNLOAD VARIABLE FILES

ds_grid = load_grid(datdir)

year = startyr
monthstrs, yearstrs = [], []

i = get_starting_i(startmo)

#Iterate over all specified months
while i <= mos:
    
    monthstr, yearstr = get_monthstr(i), str(year)
    endmonth = get_month_end(monthstr, yearstr)
    
    monthstrs.append(monthstr)
    yearstrs.append(yearstr)
    
    StartDate, EndDate = yearstr + "-" + monthstr + "-02", yearstr + "-" + monthstr + "-" + endmonth
    
    #Download monthly-averaged velocity file
    ecco_podaac_download(ShortName=vel_monthly_shortname, StartDate=StartDate, EndDate=EndDate, download_root_dir=datdir, n_workers=6, force_redownload=False)
    
    #Download monthly-averaged density/pressure file
    ecco_podaac_download(ShortName=denspress_monthly_shortname, StartDate=StartDate, EndDate=EndDate, download_root_dir=datdir, n_workers=6, force_redownload=False)
    
    if i % 12 == 0:
        year += 1 #Go to next year
        
    i += 1 #Go to next month
        
##############################

#GET FILE LISTS

vel_dir = join(datdir, vel_monthly_shortname)
denspress_dir = join(datdir, denspress_monthly_shortname)

##############################

#CREATE MONTHLY PLOTS OF GEOSTROPHIC VELOCITY

#Iterate over all specified depths
for k in range(kmin, kmax + 1):
    
    ds_vels, ds_denspressures = [], []
    
    for m in range(mos):
        
        monthstr, yearstr = monthstrs[m], yearstrs[m]
        
        curr_vel_file = join(vel_dir, vel_monthly_nc_str+yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
        curr_denspress_file = join(denspress_dir, denspress_monthly_nc_str+yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
        
        ds_vel_mo = load_dataset(curr_vel_file) #Load monthly velocity file into workspace
        
        ds_vels.append(ds_vel_mo)
        
        #Interpolate velocities to centres of grid cells
        (ds_vel_mo['UVEL']).data, (ds_vel_mo['VVEL']).data = (ds_vel_mo['UVEL']).values, (ds_vel_mo['VVEL']).values
        
        #Load monthly density-/pressure-anomaly file into workspace
        ds_denspress_mo = load_dataset(curr_denspress_file)
        
        ds_denspressures.append(ds_denspress_mo)
        
        densanom = ds_denspress_mo.RHOAnoma #Get density data
        dens = densanom + rho_ref #Compute absolute density
        
        #Update density DataArray
        
        dens.name = 'RHO'
        dens.attrs.update({'long_name': 'In-situ seawater density', 'units': 'kg m-3'})
        
        pressanom = ds_denspress_mo.PHIHYDcR #Get pressure data
        press = rho_ref * pressanom #Quantity to differentiate
        
        #Compute geostrophic velocity
        u_g, v_g = comp_geos_vel(ds_grid, press, dens, ds_vel_mo)
    
        Delta_u = comp_delta_u_norm(ds_grid, 1, ds_vel_mo, u_g, v_g) #Compute Delta-u
         
        new_grid_lon_centers, new_grid_lat_centers, new_grid_lon_edges, new_grid_lat_edges, \
    field_nearest = ecco.resample_to_latlon(ds_grid.XC, ds_grid.YC, Delta_u.values.squeeze(), 70, 89, 0.25, \
                                            -180, -90, 0.25, fill_value=np.NaN, \
                                            mapping_method='nearest_neighbor', radius_of_influence=120000)
            
            
        fig = plt.figure()
        ax = fig.add_subplot(projection=ccrs.NorthPolarStereo())
        cs1 = ax.contourf(new_grid_lon_centers, new_grid_lat_centers, field_nearest, \
                      transform=ccrs.PlateCarree(), extend='both', cmap='viridis_r', levels=np.linspace(0.01, 1, 20))
        cbar = fig.colorbar(cs1)
        ax.set_title(r'$|\Delta u|$')
        ax.add_feature(cfeature.LAND)
        ax.coastlines()
        ax.gridlines()
        plt.savefig('test.png')