"""
Master script for Beaufort Sim visualization.

Rosalie Cormier, 2023
"""

##############################

#IMPORTS AND SETTINGS

import sys
import os
import argparse
import ecco_v4_py as ecco
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

from os.path import expanduser, join

from ecco_general import load_grid, get_monthstr, load_dataset, ds_to_field, comp_residuals, rotate_vector, get_vector_partner
from ecco_visualization import cbar_label, contourf_quiver_title, ArcCir_contourf_quiver, ArcCir_contourf_quiver_grid
from ecco_field_variables_new import get_field_vars, get_variable_str

vir_nanmasked = plt.get_cmap('viridis_r').copy()
vir_nanmasked.set_bad('black')

##############################

#PARSE COMMAND-LINE INPUT AND SET GLOBAL VARIABLES

parser = argparse.ArgumentParser(description="Plot various attributes in Beaufort Gyre",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#Spatial bounds

parser.add_argument("--lats", type=float, help="Bounding latitudes", nargs=2, \
                    default=[70.0, 85.0])
parser.add_argument("--lons", type=float, help="Bounding longitudes", nargs=2, \
                    default=[-180.0, -90.0])
parser.add_argument("--res", type=float, help="Lat/lon resolution in degrees", nargs=1, \
                    default=1.0)
parser.add_argument("--kvals", type=int, help="Bounding k-values", nargs=2, default=[0, 1])

#Temporal bounds

parser.add_argument("start", type=int, help="Start year") #This argument is required
parser.add_argument("--years", type=int, help="Total number of years to plot", default=1)
parser.add_argument("--seasonal", type=bool, help="Whether to plot specific seasons", \
                    default=False)

#Attributes

parser.add_argument("--scalar", type=str, help="Name of scalar attribute", nargs="+", default="PHIHYDcR")
parser.add_argument("--scalarECCO", type=bool, help="Whether scalar field comes from ECCO files", \
                    default=True)
parser.add_argument("--vminmax", type=float, help="Minimum/maximum scalar values", nargs=2, \
                    default=[-1, 1])
parser.add_argument("--xvec", type=str, help="Name of vector attribute (x-comp)", default=None)
parser.add_argument("--vectorECCO", type=bool, help="Whether vector field comes from ECCO files", \
                    default=True)

#Directories

parser.add_argument("--datdir", type=str, help="Directory (rel. to home) with raw ECCO data", \
                    default="Downloads")
parser.add_argument("--outdir", type=str, help="Output directory (rel. to here)", \
                    default="visualization")

args = parser.parse_args()
config = vars(args)

#Spatial bounds

latmin, latmax = config['lats'][0], config['lats'][1]
lonmin, lonmax = config['lons'][0], config['lons'][1]
resolution = config['res']
kmin, kmax = config['kvals'][0], config['kvals'][1]

#Temporal bounds

startyr = config['start']
years = config['years']
seasonal = config['seasonal']

#Attributes

scalar_attr = config['scalar']
vmin, vmax = config['vminmax'][0], config['vminmax'][1]

include_vector_field = False
xvec_attr = config['xvec']

if xvec_attr is not None:
    
    include_vector_field = True
    yvec_attr = get_vector_partner(xvec_attr)
    
    variables_str = get_variable_str(xvec_attr+yvec_attr) + "_" + get_variable_str(scalar_attr)
    
elif xvec_attr is None:
    
    variables_str = get_variable_str(scalar_attr)
    
scalarECCO, vectorECCO = config['scalarECCO'], config['vectorECCO']
    
#Directories

user_home_dir = expanduser('~')
sys.path.append(join(user_home_dir, 'ECCOv4-py'))
datdir = join(user_home_dir, config['datdir'], 'ECCO_V4r4_PODAAC')
  
outdir = join(".", config['outdir'])

if not os.path.exists(outdir):
    os.makedirs(outdir)
    
if seasonal:
    subdirs = ["seasonal"]
    
elif not seasonal:
    subdirs = ["monthly", "yearly"]

for subdir in subdirs:
    if not os.path.exists(join(outdir, subdir)):
        os.makedirs(join(outdir, subdir))
        
##############################

#GET FILE NAMES

if scalarECCO:
    
    scalar_monthly_shortname, scalar_monthly_nc_str = get_field_vars(scalar_attr)
    scalar_dir = join(datdir, scalar_monthly_shortname)
    
if include_vector_field:
    
    if vectorECCO:
    
        vector_monthly_shortname, vector_monthly_nc_str = get_field_vars(xvec_attr+yvec_attr)
        vector_dir = join(datdir, vector_monthly_shortname)

ds_grid = load_grid(datdir) #Load grid  