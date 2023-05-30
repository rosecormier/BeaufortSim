"""
Compute and plot vorticity and Okubo-Weiss parameter.

Rosalie Cormier, 2023
"""

#IMPORTS AND SETTINGS

import os
import sys
import argparse
import ecco_v4_py as ecco

from os.path import expanduser, join

from ecco_general import load_grid, get_starting_i, get_monthstr, load_dataset, get_vector_in_xy, ds_to_field, comp_temp_mean
from ecco_field_variables import get_vector_field_vars
from ecco_visualization import ArcCir_pcolormesh
from vorticity_functions import comp_vorticity

##############################

#PARSE COMMAND-LINE INPUT AND SET GLOBAL VARIABLES

parser = argparse.ArgumentParser(description="Plot vorticity and Okubo-Weiss in Beaufort Gyre", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--lats", type=float, help="Bounding latitudes", nargs=2, \
                    default=[70.0, 85.0])
parser.add_argument("--lons", type=float, help="Bounding longitudes", nargs=2, \
                    default=[-180.0, -90.0])
parser.add_argument("--month", type=str, help="Start month", default="01")
parser.add_argument("--months", type=int, help="Total number of months", default=12)
parser.add_argument("--kvals", type=int, help="Bounding k-values", nargs=2, default=[12, 13])
parser.add_argument("--res", type=float, help="Lat/lon resolution in degrees", nargs=1, default=0.25)
parser.add_argument("--datdir", type=str, help="Directory (rel. to home) to store ECCO data", default="Downloads")
parser.add_argument("--outdir", type=str, help="Output directory (rel. to here)", default="visualization")

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

outdir = join(".", config['outdir'], 'vorticity')

if not os.path.exists(outdir):
    os.makedirs(outdir)
    
#Get variables associated with velocity
vel_monthly_shortname, vel_monthly_nc_str, vel_variable = get_vector_field_vars('UVEL', 'VVEL')
    
##############################

#LOAD GRID AND GET LISTS OF MONTHS/YEARS

ds_grid = load_grid(datdir)

year = startyr
monthstrs, yearstrs = [], []

i = get_starting_i(startmo)

#Iterate over all specified months
while i <= mos:
    
    monthstr, yearstr = get_monthstr(i), str(year)
    monthstrs.append(monthstr)
    yearstrs.append(yearstr)

    if (i + 1) % 12 == 0 and (i + 1) != get_starting_i(startmo) + mos:
        year += 1 #Go to next year
        
    i += 1 #Go to next month

##############################

#GET FILE LIST
vel_dir = join(datdir, vel_monthly_shortname)

##############################

#CREATE PLOTS

#Iterate over all specified depths
for k in range(kmin, kmax + 1):
    
    UVELs, VVELs, zetas = [], [], []
    
    for m in range(mos):
        
        monthstr, yearstr = monthstrs[m], yearstrs[m]
        
        curr_vel_file = join(vel_dir, vel_monthly_nc_str+yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
        
        ds_vel_mo = load_dataset(curr_vel_file) #Load monthly velocity file into workspace
        
        #Interpolate velocities to centres of grid cells
        (ds_vel_mo['UVEL']).data, (ds_vel_mo['VVEL']).data = (ds_vel_mo['UVEL']).values, (ds_vel_mo['VVEL']).values
        UVELs.append(ds_vel_mo.UVEL)
        VVELs.append(ds_vel_mo.VVEL)
        
        xgcm_grid = ecco.get_llc_grid(ds_grid)
        
        #Compute vorticity
        ds_vel_mo['zeta'] = comp_vorticity(xgcm_grid, ds_vel_mo.UVEL, ds_vel_mo.VVEL, ds_grid.dxC, ds_grid.dyC, ds_grid.rAz, k)
        
        lon_centers, lat_centers, lon_edges, lat_edges, zeta = ds_to_field(ds_grid, ds_vel_mo, 'zeta', k, latmin, latmax, lonmin, lonmax, resolution)
        zetas.append(zeta)

    #Plot annual average vorticity
    ArcCir_pcolormesh(ds_grid, k, zetas, resolution, 'PuOr', lon_centers, lat_centers, lon_edges, lat_edges, yearstr, scalar_attr='zeta', outfile="test.png", lats_lons=[70.0, 85.0, -180.0, -90.0])
        