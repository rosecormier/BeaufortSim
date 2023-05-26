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

from ecco_general import load_grid, get_monthstr, get_month_end, load_dataset, comp_residuals, get_starting_i
from ecco_visualization import cbar_label, contourf_quiver_title, ArcCir_contourf_quiver, ArcCir_contourf_quiver_grid
from ecco_field_variables import get_scalar_field_vars, get_vector_field_vars

vir_nanmasked = plt.get_cmap('viridis_r').copy()
vir_nanmasked.set_bad('black')

##############################

#PARSE COMMAND-LINE INPUT AND SET GLOBAL VARIABLES

parser = argparse.ArgumentParser(description="Plot scalar and vector fields in Beaufort Gyre",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--lats", type=float, help="Bounding latitudes", nargs=2, \
                    default=[70.0, 85.0])
parser.add_argument("--lons", type=float, help="Bounding longitudes", nargs=2, \
                    default=[-180.0, -90.0])
parser.add_argument("--month", type=str, help="Start month", default="01")
parser.add_argument("--months", type=int, help="Total number of months", default=12)
parser.add_argument("--kvals", type=int, help="Bounding k-values", nargs=2, default=[0, 4])
parser.add_argument("--res", type=float, help="Lat/lon resolution in degrees", nargs=1, \
                    default=1.0)
parser.add_argument("--datdir", type=str, help="Directory (rel. to home) to store ECCO data", \
                    default="Downloads")
parser.add_argument("--outdir", type=str, help="Output directory (rel. to here)", \
                    default="visualization")

parser.add_argument("start", type=int, help="Start year")

args = parser.parse_args()
config = vars(args)

latmin, latmax = config['lats'][0], config['lats'][1]
lonmin, lonmax = config['lons'][0], config['lons'][1]
lats_lons = [latmin, latmax, lonmin, lonmax]
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

#Set parameters and get associated variables

scalar_attr = 'PHIHYDcR'
xvec_attr, yvec_attr = 'UVEL', 'VVEL'

vector_monthly_shortname, vector_monthly_nc_str, vector_variable = \
    get_vector_field_vars(xvec_attr, yvec_attr)
scalar_monthly_shortname, scalar_monthly_nc_str, scalar_variable = \
    get_scalar_field_vars(scalar_attr)
variables_str = vector_variable + '_' + scalar_variable
    
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
    
    #Download the monthly-averaged vector file
    ecco_podaac_download(ShortName=vector_monthly_shortname, \
                         StartDate=StartDate, EndDate=EndDate, \
                         download_root_dir=datdir, n_workers=6, 
                         force_redownload=False)
     
    #Download the monthly-averaged scalar file
    ecco_podaac_download(ShortName=scalar_monthly_shortname, \
                         StartDate=StartDate, EndDate=EndDate, \
                         download_root_dir=datdir, n_workers=6, 
                         force_redownload=False)
    
    if i % 12 == 0:
        year += 1 #Go to next year
        
    i += 1 #Go to next month

##############################

#GET FILE LISTS

vector_dir = join(datdir, vector_monthly_shortname)
scalar_dir = join(datdir, scalar_monthly_shortname)

##############################

#CREATE MONTHLY PLOTS OF VELOCITY AND PRESSURE ANOMALY

#Iterate over all specified depths
for k in range(kmin, kmax + 1):

    ds_vectors, ds_scalars = [], []

    for m in range(mos):

        monthstr, yearstr = monthstrs[m], yearstrs[m]

        curr_vector_file = join(vector_dir, vector_monthly_nc_str+yearstr+ \
                                "-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
        curr_scalar_file = join(scalar_dir, scalar_monthly_nc_str+yearstr+ \
                                "-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")

        #Load monthly vector file into workspace
        ds_vector_mo = load_dataset(curr_vector_file)

        ds_vectors.append(ds_vector_mo) 

        #Interpolate vectors to centres of grid cells
        (ds_vector_mo[xvec_attr]).data, (ds_vector_mo[yvec_attr]).data = \
            (ds_vector_mo[xvec_attr]).values, (ds_vector_mo[yvec_attr]).values

        #Load monthly scalar file into workspace
        ds_scalar_mo = load_dataset(curr_scalar_file)

        ds_scalars.append(ds_scalar_mo)

        #Plot vector and scalar fields
        ArcCir_contourf_quiver(ds_grid, k, [ds_scalar_mo], [ds_vector_mo], \
                               scalar_attr, xvec_attr, yvec_attr, resolution, \
                               vir_nanmasked, [93, 97], yearstr+"-"+monthstr, \
                               outfile=join(outdir, \
                                            '{}_k{}_{}-{}.pdf'.format(variables_str, \
                                                                      str(k), \
                                                                      monthstr, \
                                                                      yearstr)), \
                               lats_lons=lats_lons)

    #Plot all months
    ArcCir_contourf_quiver_grid(ds_grid, k, ds_scalars, ds_vectors, scalar_attr, \
                                [93, 97], xvec_attr, yvec_attr, resolution, vir_nanmasked, \
                                monthstrs, yearstrs, \
                                outfile=join(outdir, '{}_k{}_all{}.png'.format(variables_str, \
                                                                               str(k), \
                                                                               yearstr)), \
                                lats_lons=lats_lons)

    #Plot annual averages
    scalar_mean, vector_mean = ArcCir_contourf_quiver(ds_grid, k, ds_scalars, ds_vectors, \
                                                      scalar_attr, xvec_attr, yvec_attr, \
                                                      resolution, vir_nanmasked, [93, 97], \
                                                      yearstrs[0]+" average", \
                                                      outfile=join(outdir, \
                                                                   '{}_k{}_avg{}.pdf'.format(variables_str, \
                                                                                             str(k), \
                                                                                             yearstr)), \
                                                      lats_lons=lats_lons)

    #Compute residuals of monthly averages

    scalar_residuals = comp_residuals([scalar_attr], ds_scalars, scalar_mean)
    vector_residuals = comp_residuals([xvec_attr, yvec_attr], ds_vectors, vector_mean)

    #Plot residuals for all months
    ArcCir_contourf_quiver_grid(ds_grid, k, scalar_residuals, vector_residuals, \
                                scalar_attr, [-2, 2], xvec_attr, yvec_attr, resolution, \
                                'seismic', monthstrs, yearstrs, \
                                outfile=join(outdir, '{}_k{}_resids_all{}.pdf'.format(variables_str, \
                                                                                      str(k), \
                                                                                      yearstr)), \
                                lats_lons=lats_lons, resid=True)