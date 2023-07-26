"""
ECCO download routine.

Rosalie Cormier, 2023
"""

import os
import sys
import argparse

from os.path import expanduser, join
from ecco_download import ecco_podaac_download

from functions_ecco_general import get_monthstr, get_month_end, get_month_name, get_vector_partner
from functions_field_variables import get_field_vars, get_variable_str

##############################

def main(**kwargs):
    
    if not kwargs:
    
        parser = argparse.ArgumentParser(description="Download new ECCO data files", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        parser.add_argument("--month", type=str, help="Start month", default="01")
        parser.add_argument("--months", type=int, help="Total number of months", default=12)
        parser.add_argument("--variable_strs", type=str, help="Variables to download", default=[])
        parser.add_argument("--datdir", type=str, help="Directory (rel. to home) to store ECCO data", default="Downloads")

        parser.add_argument("start", type=int, help="Start year")
    
        args = parser.parse_args()
        config = vars(args)

        startmo, startyr, mos = config['month'], config['start'], config['months']
        scalars, xvectors = config['scalars'], config['xvectors']
        
        datdirshort = config['datdir']
        
    elif kwargs:
        
        startmo = kwargs.get('startmo')
        startyr = kwargs.get('startyr')
        mos = kwargs.get('months')
        datdirshort = kwargs.get('datdir')
        variable_strs = kwargs.get('variable_strs')

    homedir = expanduser('~')
    datdir = join(homedir, datdirshort, 'ECCO_V4r4_PODAAC')

    if not os.path.exists(datdir):
        os.makedirs(datdir)

    ##############################

    #DOWNLOAD VARIABLE FILES

    year = startyr

    start_i, monthname = get_month_name(startmo)
    i = start_i
    
    while i <= start_i + mos: #Iterate over all specified months

        monthstr, yearstr = get_monthstr(i), str(year)
        endmonth = get_month_end(monthstr, yearstr)

        StartDate, EndDate = yearstr + "-" + monthstr + "-02", yearstr + "-" + monthstr + "-" + endmonth

        for variable_str in variable_strs:
            
            monthly_shortname, monthly_nc_str = get_field_vars(variable_str)
            ecco_podaac_download(ShortName=monthly_shortname, StartDate=StartDate, EndDate=EndDate, download_root_dir=datdir, n_workers=6, force_redownload=False)

        if (i + 1) % 12 == 0 and (i + 1) != start_i + mos:
            year += 1 #Go to next year

        i += 1 #Go to next month
        
##############################

if __name__ == "__main__":
    main()