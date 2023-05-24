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
from ecco_field_variables import *

##############################

#PARSE COMMAND-LINE INPUT AND SET GLOBAL VARIABLES

parser = argparse.ArgumentParser(description="Plot scalar and vector fields in Beaufort Gyre",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--lats", type=float, help="Bounding latitudes", nargs=2, default=[70.0, 85.0])
parser.add_argument("--lons", type=float, help="Bounding longitudes", nargs=2, default=[-180.0, -90.0])
parser.add_argument("--month", type=str, help="Start month", default="01")
parser.add_argument("--months", type=int, help="Total number of months", default=12)
parser.add_argument("--kvals", type=int, help="Bounding k-values", nargs=2, default=[0, 4])
parser.add_argument("--res", type=float, help="Lat/lon resolution in degrees", nargs=1, default=1.0)
parser.add_argument("--datdir", type=str, help="Directory (rel. to home) to store ECCO data", default="Downloads")
parser.add_argument("--outdir", type=str, help="Output directory (rel. to here)", default="visualization")

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
    
month_dict = {0: "01", 1: "02", 2: "03", 3: "04", 4: "05", 5: "06",
              6: "07", 7: "08", 8: "09", 9: "10", 10: "11", 11: "12"}
month_key_list, month_val_list = list(month_dict.keys()), list(month_dict.values())

#Set parameters and get associated variables

k = 1
scalar_attr = 'PHIHYDcR'
xvec_attr = 'UVEL'
yvec_attr = 'VVEL'

vector_monthly_shortname, vector_monthly_nc_str, vector_variable = get_vector_field_vars(xvec_attr, yvec_attr)
scalar_monthly_shortname, scalar_monthly_nc_str, scalar_variable = get_scalar_field_vars(scalar_attr)
variables_str = vector_variable + '_' + scalar_variable
    
##############################

#DOWNLOAD ECCO FILES FOR SPECIFIED MONTHS

grid_params_shortname = "ECCO_L4_GEOMETRY_LLC0090GRID_V4R4"

month_index = month_val_list.index(startmo)
i = month_key_list[month_index]
year = startyr
monthstrs, yearstrs = [], []

#Iterate over all specified months
while i < mos:

    monthstr = month_dict[i % 12]
    yearstr = str(year)
    endmonth = get_month_end(monthstr, yearstr)
    
    monthstrs.append(monthstr)
    yearstrs.append(yearstr)
    
    #Download the monthly-averaged velocity file
    ecco_podaac_download(ShortName=vector_monthly_shortname, StartDate=yearstr+"-"+monthstr+"-02", 
                         EndDate=yearstr+"-"+monthstr+"-"+endmonth, download_root_dir=datdir, n_workers=6, 
                         force_redownload=False)
     
    #Download the monthly-averaged density-/pressure-anomaly file
    ecco_podaac_download(ShortName=scalar_monthly_shortname, StartDate=yearstr+"-"+monthstr+"-02", 
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

vector_dir = join(datdir, vector_monthly_shortname)
scalar_dir = join(datdir, scalar_monthly_shortname)

grid_params_file = "GRID_GEOMETRY_ECCO_V4r4_native_llc0090.nc"
grid_params_file_path = join(datdir, grid_params_shortname, grid_params_file)

#Load grid parameters
ds_grid = xr.open_dataset(grid_params_file_path)

##############################

#CREATE MONTHLY PLOTS OF VELOCITY AND PRESSURE ANOMALY

vir_nanmasked = plt.get_cmap('viridis_r').copy()
vir_nanmasked.set_bad('black')

ds_vectors, ds_scalars = [], []

for m in range(mos):
    
    monthstr = monthstrs[m]
    yearstr = yearstrs[m]
    
    curr_vector_file = join(vector_dir, vector_monthly_nc_str+yearstr+
                         "-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
    curr_scalar_file = join(scalar_dir, scalar_monthly_nc_str+
                               yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")

    #Load monthly velocity file into workspace
    ds_vector_mo = xr.open_mfdataset(curr_vector_file, parallel=True, data_vars='minimal', coords='minimal', 
                              compat='override')
    
    ds_vectors.append(ds_vector_mo) 
    
    #Interpolate velocities to centres of grid cells
    (ds_vector_mo[xvec_attr]).data, (ds_vector_mo[yvec_attr]).data = (ds_vector_mo[xvec_attr]).values, (ds_vector_mo[yvec_attr]).values
    
    #Load monthly density-/pressure-anomaly file into workspace
    ds_scalar_mo = xr.open_mfdataset(curr_scalar_file, parallel=True, data_vars='minimal', coords='minimal', 
                                    compat='override')
    
    ds_scalars.append(ds_scalar_mo)

    #Plot vector and scalar fields
    ArcCir_contourf_quiver(ds_grid, k, [ds_scalar_mo], [ds_vector_mo], scalar_attr, xvec_attr, yvec_attr, resolution, vir_nanmasked, [93, 97], yearstr+"-"+monthstr, outfile=join(outdir, '{}_{}-{}.pdf'.format(variables_str, monthstr, yearstr)), lats_lons=lats_lons)

#Plot all months
ArcCir_contourf_quiver_grid(ds_grid, k, ds_scalars, ds_vectors, scalar_attr, [93, 97], xvec_attr, yvec_attr, resolution, 
                           vir_nanmasked, monthstrs, yearstrs, outfile=join(outdir, 'u_p_anom_all{}.png'.format(yearstr)),
                           lats_lons=lats_lons)

#Plot annual averages
scalar_mean, vector_mean = ArcCir_contourf_quiver(ds_grid, k, ds_scalars, ds_vectors, scalar_attr, xvec_attr, yvec_attr, resolution, vir_nanmasked, [93, 97], yearstrs[0]+" average", outfile=join(outdir, '{}_avg{}.pdf'.format(variables_str, yearstr)), lats_lons=lats_lons)

#Compute residuals of monthly averages

scalar_residuals = comp_residuals([scalar_attr], ds_scalars, scalar_mean)
vector_residuals = comp_residuals([xvec_attr, yvec_attr], ds_vectors, vector_mean)

#Plot residuals for all months
ArcCir_contourf_quiver_grid(ds_grid, k, scalar_residuals, vector_residuals, scalar_attr, [-2, 2], xvec_attr, yvec_attr, resolution, 'seismic', monthstrs, yearstrs, outfile=join(outdir, '{}_resids_all{}.pdf'.format(variables_str, yearstr)), lats_lons=lats_lons)