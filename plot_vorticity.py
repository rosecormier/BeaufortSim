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

from ecco_general import load_grid, get_starting_i, get_monthstr, load_dataset, get_vector_in_xy
from ecco_field_variables import get_vector_field_vars
from vorticity_functions import comp_vorticity

##############################

#PARSE COMMAND-LINE INPUT AND SET GLOBAL VARIABLES

parser = argparse.ArgumentParser(description="Plot vorticity and Okubo-Weiss in Beaufort Gyre", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--month", type=str, help="Start month", default="01")
parser.add_argument("--months", type=int, help="Total number of months", default=12)
parser.add_argument("--kvals", type=int, help="Bounding k-values", nargs=2, default=[12, 13])
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
    
    for m in range(mos):
        
        monthstr, yearstr = monthstrs[m], yearstrs[m]
        
        curr_vel_file = join(vel_dir, vel_monthly_nc_str+yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
        
        ds_vel_mo = load_dataset(curr_vel_file) #Load monthly velocity file into workspace
        
        #Interpolate velocities to centres of grid cells
        
        (ds_vel_mo['UVEL']).data, (ds_vel_mo['VVEL']).data = (ds_vel_mo['UVEL']).values, (ds_vel_mo['VVEL']).values
        velocity_interp = get_vector_in_xy(ds_grid, ds_vel_mo, 'UVEL', 'VVEL') 
        u, v = velocity_interp['X'], velocity_interp['Y']
        #u, v = (u.isel(k=k)).squeeze(), (v.isel(k=k)).squeeze()

        xgcm_grid = ecco.get_llc_grid(ds_grid)
        
        print(comp_vorticity(xgcm_grid, ds_vel_mo.UVEL, ds_vel_mo.VVEL, ds_grid.dxC, ds_grid.dyC, ds_grid.rAz, k))