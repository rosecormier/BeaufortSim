"""
Rosalie Cormier, 2023, based on code by Andrew Delman
"""

##############################

#IMPORTS AND SETTINGS

import os
import sys
import argparse

from os.path import expanduser, join

from ecco_field_variables import *

##############################

#PARSE COMMAND-LINE INPUT AND SET GLOBAL VARIABLES

parser = argparse.ArgumentParser(description="Plot geostrophic velocity in Beaufort Gyre", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--datdir", type=str, help="Directory (rel. to home) to store ECCO data", default="Downloads")
parser.add_argument("--outdir", type=str, help="Output directory (rel. to here)", default="visualization")

args = parser.parse_args()
config = vars(args)

user_home_dir = expanduser('~')
sys.path.append(join(user_home_dir, 'ECCOv4-py'))
datdir = join(user_home_dir, config['datdir'], 'ECCO_V4r4_PODAAC')

if not os.path.exists(datdir):
    os.makedirs(datdir)
    
outdir = join(".", config['outdir'], 'geostrophic')

if not os.path.exists(outdir):
    os.makedirs(outdir)
    
#Get variables associated with velocity
vel_monthly_shortname, vel_monthly_nc_str, vel_variable = get_vector_field_vars('UVEL', 'VVEL')

##############################

#LOAD GRID AND DOWNLOAD VARIABLE FILES

ds_grid = load_grid(datdir)

month_index = 