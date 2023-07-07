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

from ecco_general import load_grid, get_monthstr, load_dataset, ds_to_field, comp_residuals, rotate_vector, get_vector_partner, ecco_resample, get_season_months_and_years, get_scalar_in_xy, get_vector_in_xy
from ecco_visualization import ArcCir_pcolormesh, ArcCir_contourf_quiver, ArcCir_contourf_quiver_grid
from ecco_field_variables import get_field_vars, get_variable_str
from geostrophic_functions import rotate_u_g, comp_geos_metric
from vorticity_functions import comp_local_Ro, comp_normal_strain, comp_shear_strain, comp_OkuboWeiss

#The following are scripts that are imported as modules but may be run within this script

import compute_monthly_avgs
import download_new_data
import save_annual_avgs
import save_seasonal_avgs

##############################

def get_parser():

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
    parser.add_argument("--seasonal", dest='seasonal', help="Whether to plot specific seasons", \
                            default=False, action='store_true')
    parser.add_argument("--seasonmonths", type=str, help="Start and end months of season", nargs=2, \
                            default=["01", "01"])

    #Attributes

    parser.add_argument("--scalar", type=str, help="Name of scalar attribute", default="ZETA")
    parser.add_argument("--scalarECCO", dest='scalarECCO', help="Whether scalar field comes from ECCO files", \
                            default=False, action='store_true')
    parser.add_argument("--vminmax", type=float, help="Minimum/maximum scalar values", nargs=2, \
                            default=[1, 1])
    parser.add_argument("--xvec", type=str, help="Name of vector attribute (x-comp)", default=None)
    parser.add_argument("--vectorECCO", dest='vectorECCO', help="Whether vector field comes from ECCO files", \
                            default=False, action='store_true')

    #Directories

    parser.add_argument("--datdir", type=str, help="Directory (rel. to home) with raw ECCO data", \
                        default="Downloads")
    parser.add_argument("--compdatdir", type=str, help="Directory (rel. to here) with computed monthly data", \
                        default="computed_monthly")
    parser.add_argument("--seasonaldatdir", type=str, help="Directory (rel. to here) with seasonal avgs", \
                        default="seasonal_averages")
    parser.add_argument("--yearlydatdir", type=str, help="Directory (rel. to here) with annual avgs", \
                        default="yearly_averages")
    parser.add_argument("--outdir", type=str, help="Output directory (rel. to here)", \
                        default="visualization")

    #Visualization
    parser.add_argument("--cmap", type=str, help="MPL colormap", default="viridis_r")
    
    return parser

def scalarECCO_load_dataset(scalar_dir, scalar_monthly_nc_str, yearstr, monthstr, year, scalar_attr, datdir):
    
    curr_scalar_file = join(scalar_dir, scalar_monthly_nc_str+yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")

    if not os.path.exists(curr_scalar_file):
                                
        #If the file doesn't exist, download it
        download_new_data.main(startmo="01", startyr=year, months=12, scalars=[scalar_attr], xvectors=None, datdir=datdir)
                            
    ds_scalar_mo = load_dataset(curr_scalar_file) #Load monthly scalar file into workspace
    
    return ds_scalar_mo

def main():

    vir_nanmasked = plt.get_cmap('viridis_r').copy()
    vir_nanmasked.set_bad('black')

    ##############################
    
    parser = get_parser()
    args = parser.parse_args()
    config = vars(args)

    #Spatial bounds

    latmin, latmax = config['lats'][0], config['lats'][1]
    lonmin, lonmax = config['lons'][0], config['lons'][1]
    lats_lons = [latmin, latmax, lonmin, lonmax]

    resolution = config['res']
    kmin, kmax = config['kvals'][0], config['kvals'][1]

    #Temporal bounds

    startyr = config['start']
    years = config['years']
    seasonal = config['seasonal']

    if seasonal:
        season_start, season_end = config['seasonmonths']

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

    homedir = expanduser('~')
    sys.path.append(join(homedir, 'ECCOv4-py'))
    datdir = join(homedir, config['datdir'], 'ECCO_V4r4_PODAAC')

    compdatdir = join(".", config['compdatdir'])

    outdir = join(".", config['outdir'])

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if seasonal:

        subdirs = ["seasonal", "interannual"]

        seasonaldatdir = join(".", config['seasonaldatdir'])

    elif not seasonal:

        subdirs = ["monthly", "yearly"]

        yearlydatdir = join(".", config['yearlydatdir'])

    for subdir in subdirs:
        if not os.path.exists(join(outdir, subdir)):
            os.makedirs(join(outdir, subdir))

    cmap = config['cmap']

    ##############################

    #GET FILE NAMES

    scalar_monthly_shortname, scalar_monthly_nc_str = get_field_vars(scalar_attr)

    if scalarECCO:
        scalar_dir = join(datdir, scalar_monthly_shortname)

    elif not scalarECCO:
        scalar_dir = join(compdatdir, scalar_monthly_shortname)

    if include_vector_field:

        vector_monthly_shortname, vector_monthly_nc_str = get_field_vars(xvec_attr+yvec_attr)

        if vectorECCO:
            vector_dir = join(datdir, vector_monthly_shortname)

        elif not vectorECCO:
            vector_dir = join(compdatdir, vector_monthly_shortname)

    ds_grid = load_grid(datdir) #Load grid  

    ##############################

    #CASE WHERE ONLY A SCALAR IS PROVIDED

    if not include_vector_field:

        for k in range(kmin, kmax + 1): #Iterate over specified depths

            if not seasonal: #Case where we plot every month

                for i in range(years): #Iterate over specified years

                    year = startyr + i
                    yearstr = str(year)
                    
                    if scalar_attr == 'ZETA':
                        Ro_l_list, OW_list = [], []

                    for m in range(12): #Iterate over months
                        
                        monthstr = get_monthstr(m)

                        if scalarECCO:
                            ds_scalar_mo = scalarECCO_load_dataset(scalar_dir, scalar_monthly_nc_str, yearstr, monthstr, year, scalar_attr, datdir=config['datdir'])

                        elif not scalarECCO:

                            curr_scalar_file = join(scalar_dir, scalar_monthly_nc_str+yearstr+"-"+monthstr+".nc")

                            if os.path.exists(curr_scalar_file): #Look for the file
                                ds_scalar_mo = xr.open_mfdataset(curr_scalar_file, engine="scipy") #Load monthly scalar file into workspace

                            else: #If it doesn't exist, compute it
                                
                                compute_monthly_avgs.main(latmin=70.0, latmax=85.0, lonmin=-180.0, lonmax=-90.0, startyr=year, years=1, datdir=config['datdir'], outdir=compdatdir)
                                ds_scalar_mo = load_dataset(curr_scalar_file)
                                
                            ds_grid = get_scalar_in_xy(ds_grid, ds_scalar_mo, scalar_attr)

                        #Convert scalar DataSet to useful field
                        lon_centers, lat_centers, lon_edges, lat_edges, scalar = ds_to_field(ds_grid, ds_scalar_mo.isel(k=k), scalar_attr, latmin, latmax, lonmin, lonmax, resolution)

                        #Plot monthly scalar data
                        ArcCir_pcolormesh(ds_grid, k, [scalar], resolution, cmap, lon_centers, lat_centers, monthstr+"-"+yearstr, scalar_attr, scalar_bounds=[vmin, vmax], extend='both', outfile=join(outdir, 'monthly', '{}_k{}_{}{}.pdf'.format(variables_str, str(k), monthstr, yearstr)), lats_lons=lats_lons) 

                        if scalar_attr == 'ZETA': #If vorticity, also compute and plot Ro_l, OW
                            
                            Ro_l = comp_local_Ro(scalar, lat_centers) #Compute Ro_l
                            
                            #Plot Ro_l
                            ArcCir_pcolormesh(ds_grid, k, [Ro_l], resolution, 'Reds', lon_centers, lat_centers, monthstr+"-"+yearstr, 'Ro_l', scalar_bounds=[0, 1], extend='max', outfile=join(outdir, 'monthly', 'localRo_k{}_{}{}.pdf'.format(str(k), monthstr, yearstr)), lats_lons=lats_lons)
                            
                            Ro_l_list.append(Ro_l)
                            
                            #Get monthly velocity data
                            vel_monthly_shortname, vel_monthly_nc_str = get_field_vars('UVELVVEL')
                                
                            ds_vel_mo = scalarECCO_load_dataset(join(datdir, vel_monthly_shortname), vel_monthly_nc_str, yearstr, monthstr, year, 'UVEL', datdir=config['datdir'])
                            (ds_vel_mo['UVEL']).data, (ds_vel_mo['VVEL']).data = (ds_vel_mo['UVEL']).values, (ds_vel_mo['VVEL']).values
                            
                            #Compute strain terms
                            
                            xgcm_grid = ecco.get_llc_grid(ds_grid)
                            normal_strain = comp_normal_strain(xgcm_grid, ds_vel_mo['UVEL'], ds_vel_mo['VVEL'], ds_grid.dxG, ds_grid.dyG, ds_grid.rA).isel(k=k).squeeze()
                            shear_strain = comp_shear_strain(xgcm_grid, ds_vel_mo['UVEL'], ds_vel_mo['VVEL'], ds_grid.dxC, ds_grid.dyC, ds_grid.rAz).isel(k=k).squeeze()
                            normal_strain= ecco_resample(ds_grid, normal_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
                            shear_strain = ecco_resample(ds_grid, shear_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
                            
                            OW = comp_OkuboWeiss(scalar, normal_strain, shear_strain) #Compute OW
                            
                            #Plot OW
                            ArcCir_pcolormesh(ds_grid, k, [OW], resolution, 'seismic', lon_centers, lat_centers, monthstr+"-"+yearstr, 'OW', scalar_bounds=[-1e-13, 1e-13], extend='both', outfile=join(outdir, 'monthly', 'OW_k{}_{}{}.pdf'.format(str(k), monthstr, yearstr)), lats_lons=lats_lons)
                            
                            OW_list.append(OW)
       
                    #Get annually-averaged data

                    scalar_annual_file = join(yearlydatdir, "avg_"+scalar_attr+"_"+yearstr+".nc")
                    
                    if not os.path.exists(scalar_annual_file): #If it doesn't exist, compute it
                        
                        if scalarECCO:
                            datdirshort, usecompdata = 'Downloads', False
                            
                        elif not scalarECCO:
                            datdirshort, usecompdata = 'computed_monthly', True
                            
                        save_annual_avgs.main(years=[year], field=scalar_attr, datdir=config['datdir'], usecompdata=usecompdata, outdir=yearlydatdir)
                    
                    ds_scalar_year = xr.open_mfdataset(scalar_annual_file, engine="scipy")
                    ds_scalar_year.load()

                    #Convert scalar DataSet to useful field
                    lon_centers, lat_centers, lon_edges, lat_edges, scalar_year = ds_to_field(ds_grid, ds_scalar_year.isel(k=k), scalar_attr, latmin, latmax, lonmin, lonmax, resolution)

                    #Plot annual average
                    ArcCir_pcolormesh(ds_grid, k, [scalar_year], resolution, cmap, lon_centers, lat_centers, yearstr, scalar_attr, scalar_bounds=[vmin, vmax], extend='both', outfile=join(outdir, 'yearly', '{}_k{}_{}.pdf'.format(variables_str, str(k), yearstr)), lats_lons=lats_lons) 
                    
                    if scalar_attr == 'ZETA': #If vorticity, also compute and plot annual Ro_l, OW
                        
                        ArcCir_pcolormesh(ds_grid, k, Ro_l_list, resolution, 'Reds', lon_centers, lat_centers, yearstr, 'Ro_l', scalar_bounds=[0, 1], extend='max', outfile=join(outdir, 'yearly', 'localRo_k{}_{}.pdf'.format(str(k), yearstr)), lats_lons=lats_lons)
                        
                        ArcCir_pcolormesh(ds_grid, k, OW_list, resolution, 'seismic', lon_centers, lat_centers, yearstr, 'OW', scalar_bounds=[-1e-13, 1e-13], extend='both', outfile=join(outdir, 'yearly', 'OW_k{}_{}.pdf'.format(str(k), yearstr)), lats_lons=lats_lons)

            elif seasonal: #Case where we plot one season per year

                season_months, season_years = get_season_months_and_years(season_start, season_end)
                seas_monthstr = season_months[0] + "-" + season_months[-1] #For titles

                data_seasons = []
                
                if scalar_attr == 'ZETA':
                        Ro_l_list, OW_list = [], []

                for i in range(years): #Iterate over specified years

                    year = startyr + i
                    yearstr = str(year)
                    endyearstr = str(year + season_years[-1])

                    #Get seasonally-averaged data

                    scalar_seas_file = join(seasonaldatdir, "avg_"+scalar_attr+"_"+season_start+yearstr+"-"+season_end+endyearstr+".nc")
                    
                    if not os.path.exists(scalar_seas_file): #If it doesn't exist, compute it
                        
                        if scalarECCO:
                            datdirshort, usecompdata = 'Downloads', False
                            
                        elif not scalarECCO:
                            datdirshort, usecompdata = 'computed_monthly', True
                        
                        save_seasonal_avgs.main(field=scalar_attr, years=[year], start_month=season_start, end_month=season_end, usecompdata=usecompdata, datdir=datdirshort, outdir=seasonaldatdir)
                    
                    ds_scalar_seas = xr.open_mfdataset(scalar_seas_file, engine="scipy")
                    ds_scalar_seas.load()

                    #Convert scalar DataSet to useful field
                    lon_centers, lat_centers, lon_edges, lat_edges, scalar_seas = ds_to_field(ds_grid, ds_scalar_seas.isel(k=k), scalar_attr, latmin, latmax, lonmin, lonmax, resolution)

                    seas_yearstr = yearstr + "-" + str(year + season_years[-1]) #For titles

                    #Plot seasonal average
                    ArcCir_pcolormesh(ds_grid, k, [scalar_seas], resolution, cmap, lon_centers, lat_centers, '{}, {}'.format(seas_monthstr, seas_yearstr), scalar_attr, scalar_bounds=[vmin, vmax], extend='both', outfile=join(outdir, 'seasonal', '{}_k{}_{}_{}.pdf'.format(variables_str, str(k), seas_monthstr, seas_yearstr)), lats_lons=lats_lons)

                    data_seasons.append(scalar_seas)
                    
                    if scalar_attr == 'ZETA': #If vorticity, also compute and plot Ro_l, OW
                            
                        Ro_l = comp_local_Ro(scalar_seas, lat_centers) #Compute Ro_l
                            
                        #Plot Ro_l
                        ArcCir_pcolormesh(ds_grid, k, [Ro_l], resolution, 'Reds', lon_centers, lat_centers, '{}, {}'.format(seas_monthstr, seas_yearstr), 'Ro_l', scalar_bounds=[0, 1], extend='max', outfile=join(outdir, 'seasonal', 'localRo_k{}_{}_{}.pdf'.format(str(k), seas_monthstr, seas_yearstr)), lats_lons=lats_lons)
                            
                        Ro_l_list.append(Ro_l)
                            
                        vel_seas_file = join(seasonaldatdir, "avg_UVELVVEL_"+season_start+yearstr+"-"+season_end+endyearstr+".nc")
                            
                        if not os.path.exists(vel_seas_file): #If it doesn't exist, compute it
                            save_seasonal_avgs.main(field='UVELVVEL', years=[year], start_month=season_start, end_month=season_end, usecompdata=False, datdir='Downloads', outdir=seasonaldatdir)
                                
                        ds_vel_seas = xr.open_mfdataset(vel_seas_file, engine="scipy")
                        ds_vel_seas.load()

                        #Compute strain terms
                            
                        xgcm_grid = ecco.get_llc_grid(ds_grid)
                        normal_strain = comp_normal_strain(xgcm_grid, ds_vel_seas['UVEL'], ds_vel_seas['VVEL'], ds_grid.dxG, ds_grid.dyG, ds_grid.rA).isel(k=k).squeeze()
                        shear_strain = comp_shear_strain(xgcm_grid, ds_vel_seas['UVEL'], ds_vel_seas['VVEL'], ds_grid.dxC, ds_grid.dyC, ds_grid.rAz).isel(k=k).squeeze()
                        normal_strain= ecco_resample(ds_grid, normal_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
                        shear_strain = ecco_resample(ds_grid, shear_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
                            
                        OW = comp_OkuboWeiss(scalar_seas, normal_strain, shear_strain) #Compute OW
                            
                        #Plot OW
                        ArcCir_pcolormesh(ds_grid, k, [OW], resolution, 'seismic', lon_centers, lat_centers, '{}, {}'.format(seas_monthstr, seas_yearstr), 'OW', scalar_bounds=[-1e-13, 1e-13], extend='both', outfile=join(outdir, 'seasonal', 'OW_k{}_{}_{}.pdf'.format(str(k), seas_monthstr, seas_yearstr)), lats_lons=lats_lons)
                            
                        OW_list.append(OW)
                
                if years != 1: #If there is more than one season to average over
                
                    seas_yearstr = str(startyr) + "-" + str(startyr + (years-1) + season_years[-1]) #For titles

                    #Plot average over all seasons
                    ArcCir_pcolormesh(ds_grid, k, data_seasons, resolution, cmap, lon_centers, lat_centers, '{}, {}'.format(seas_monthstr, seas_yearstr), scalar_attr, scalar_bounds=[vmin, vmax], extend='both', outfile=join(outdir, 'interannual', '{}_k{}_{}_{}.pdf'.format(variables_str, str(k), seas_monthstr, seas_yearstr)), lats_lons=lats_lons)

                    if scalar_attr == 'ZETA': #If vorticity, also compute and plot interannual Ro_l, OW

                        ArcCir_pcolormesh(ds_grid, k, Ro_l_list, resolution, 'Reds', lon_centers, lat_centers, '{}, {}'.format(seas_monthstr, seas_yearstr), 'Ro_l', scalar_bounds=[0, 1], extend='max', outfile=join(outdir, 'interannual', 'localRo_k{}_{}_{}.pdf'.format(str(k), seas_monthstr, seas_yearstr)), lats_lons=lats_lons)

                        ArcCir_pcolormesh(ds_grid, k, OW_list, resolution, 'seismic', lon_centers, lat_centers, '{}, {}'.format(seas_monthstr, seas_yearstr), 'OW', scalar_bounds=[-1e-13, 1e-13], extend='both', outfile=join(outdir, 'interannual', 'OW_k{}_{}_{}.pdf'.format(str(k), seas_monthstr, seas_yearstr)), lats_lons=lats_lons)

    ##############################

    #CASE WHERE VECTOR AND SCALAR ARE PROVIDED

    elif include_vector_field:

        for k in range(kmin, kmax + 1): #Iterate over specified depths

            if not seasonal: #Case where we plot every month

                for i in range(years): #Iterate over specified years

                    year = startyr + i
                    yearstr = str(year)

                    for m in range(12): #Iterate over months

                        monthstr = get_monthstr(m)

                        if scalarECCO:
                            ds_scalar_mo = scalarECCO_load_dataset(scalar_dir, scalar_monthly_nc_str, yearstr, monthstr, year, scalar_attr, datdir=config['datdir'])

                        elif not scalarECCO:

                            curr_scalar_file = join(scalar_dir, scalar_monthly_nc_str+yearstr+"-"+monthstr+".nc")

                            if os.path.exists(curr_scalar_file): #Look for the file
                                ds_scalar_mo = xr.open_mfdataset(curr_scalar_file, engine="scipy") #Load monthly scalar file into workspace

                            else: #If it doesn't exist, compute it
                                
                                compute_monthly_avgs.main(latmin=70.0, latmax=85.0, lonmin=-180.0, lonmax=-90.0, startyr=year, years=1, datdir=config['datdir'], outdir=compdatdir)
                                ds_scalar_mo = load_dataset(curr_scalar_file)
                                
                            ds_grid = get_scalar_in_xy(ds_grid, ds_scalar_mo, scalar_attr)

                        #Convert scalar DataSet to useful field
                        lon_centers, lat_centers, lon_edges, lat_edges, scalar = ds_to_field(ds_grid, ds_scalar_mo.isel(k=k), scalar_attr, latmin, latmax, lonmin, lonmax, resolution)
                        
                        if vectorECCO:
                            
                            ds_vector_mo = scalarECCO_load_dataset(vector_dir, vector_monthly_nc_str, yearstr, monthstr, year, xvec_attr, datdir=config['datdir'])
                            
                            #Interpolate and rotate vector
                            
                            (ds_vector_mo[xvec_attr]).data, (ds_vector_mo[yvec_attr]).data = (ds_vector_mo[xvec_attr]).values, (ds_vector_mo[yvec_attr]).values
                            vecE, vecN = rotate_vector(ds_grid, ds_vector_mo, xvec_attr, yvec_attr)
                            vecE, vecN = vecE.isel(k=k).squeeze(), vecN.isel(k=k).squeeze()
                            
                        elif not vectorECCO:
                            
                            curr_vector_file = join(vector_dir, vector_monthly_nc_str+yearstr+"-"+monthstr+".nc")

                            if os.path.exists(curr_vector_file): #Look for the file
                                ds_vector_mo = xr.open_mfdataset(curr_vector_file, engine="scipy") #Load monthly vector file into workspace

                            else: #If it doesn't exist, compute it
                                
                                compute_monthly_avgs.main(latmin=70.0, latmax=85.0, lonmin=-180.0, lonmax=-90.0, startyr=year, years=1, datdir=config['datdir'], outdir=compdatdir)
                                ds_vector_mo = load_dataset(curr_vector_file)

                            vecE, vecN = rotate_u_g(ds_grid, ds_vector_mo[xvec_attr], ds_vector_mo[yvec_attr], k)
                            
                            if xvec_attr == 'UG': #If u_g, also plot geostrophy metric
                                
                                vel_monthly_shortname, vel_monthly_nc_str = get_field_vars('UVELVVEL')
                                
                                ds_vel_mo = scalarECCO_load_dataset(join(datdir, vel_monthly_shortname), vel_monthly_nc_str, yearstr, monthstr, year, 'UVEL', datdir=config['datdir'])
                                
                                #Interpolate velocities to centres
                                (ds_vel_mo['UVEL']).data, (ds_vel_mo['VVEL']).data = (ds_vel_mo['UVEL']).values, (ds_vel_mo['VVEL']).values
                                velocity_interp = get_vector_in_xy(ds_grid, ds_vel_mo, 'UVEL', 'VVEL') 
                                u, v = velocity_interp['X'].isel(k=k), velocity_interp['Y'].isel(k=k)

                                #Compute our geostrophy metric
                                Delta_u = comp_geos_metric(u.squeeze(), v.squeeze(), vecE, vecN)
                                
                                lon_centers, lat_centers, lon_edges, lat_edges, Delta_u_plot = ecco_resample(ds_grid, Delta_u, latmin, latmax, lonmin, lonmax, resolution)
                                
                                ArcCir_pcolormesh(ds_grid, k, [Delta_u_plot], resolution, 'Reds', lon_centers, lat_centers, monthstr+"-"+yearstr, 'Delta_u', scalar_bounds=[0, 1], extend='max', outfile=join(outdir, 'monthly', '{}_k{}_{}{}.pdf'.format('Delta_u', str(k), monthstr, yearstr)), lats_lons=lats_lons) 
                            
                        #Plot monthly data
                        ArcCir_contourf_quiver(ds_grid, k, [scalar], [vecE], [vecN], resolution, cmap, yearstr+"-"+monthstr, lon_centers, lat_centers, scalar_attr, xvec_attr, scalar_bounds=[vmin, vmax], outfile=join(outdir, 'monthly', '{}_k{}_{}{}.pdf'.format(variables_str, str(k), monthstr, yearstr)), lats_lons=lats_lons)
                        
                        if scalar_attr == 'ZETA': #If vorticity, also compute and plot Ro_l, OW
                            
                            #Compute Ro_l
                            Ro_l = comp_local_Ro(scalar, lat_centers)
                            
                            #Plot Ro_l
                            ArcCir_pcolormesh(ds_grid, k, [Ro_l], resolution, 'Reds', lon_centers, lat_centers, monthstr+"-"+yearstr, 'Ro_l', scalar_bounds=[0, 1], extend='max', outfile=join(outdir, 'monthly', 'localRo_k{}_{}{}.pdf'.format(str(k), monthstr, yearstr)), lats_lons=lats_lons)
                            
                            #Get monthly velocity data
                            
                            vel_monthly_shortname, vel_monthly_nc_str = get_field_vars('UVELVVEL')
                                
                            ds_vel_mo = scalarECCO_load_dataset(join(datdir, vel_monthly_shortname), vel_monthly_nc_str, yearstr, monthstr, year, 'UVEL', datdir=config['datdir'])

                            (ds_vel_mo['UVEL']).data, (ds_vel_mo['VVEL']).data = (ds_vel_mo['UVEL']).values, (ds_vel_mo['VVEL']).values
                            
                            #Compute strain terms
                            
                            xgcm_grid = ecco.get_llc_grid(ds_grid)
                            normal_strain = comp_normal_strain(xgcm_grid, ds_vel_mo['UVEL'], ds_vel_mo['VVEL'], ds_grid.dxG, ds_grid.dyG, ds_grid.rA).isel(k=k).squeeze()
                            shear_strain = comp_shear_strain(xgcm_grid, ds_vel_mo['UVEL'], ds_vel_mo['VVEL'], ds_grid.dxC, ds_grid.dyC, ds_grid.rAz).isel(k=k).squeeze()
                            normal_strain= ecco_resample(ds_grid, normal_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
                            shear_strain = ecco_resample(ds_grid, shear_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
                            
                            #Compute OW
                            OW = comp_OkuboWeiss(scalar, normal_strain, shear_strain)
                            
                            #Plot OW
                            ArcCir_pcolormesh(ds_grid, k, [OW], resolution, 'seismic', lon_centers, lat_centers, monthstr+"-"+yearstr, 'OW', scalar_bounds=[-1e-13, 1e-13], extend='both', outfile=join(outdir, 'monthly', 'OW_k{}_{}{}.pdf'.format(str(k), monthstr, yearstr)), lats_lons=lats_lons)
                            
                    #Get annually-averaged data

                    scalar_annual_file = join(yearlydatdir, "avg_"+scalar_attr+"_"+yearstr+".nc")
                    vector_annual_file = join(yearlydatdir, "avg_"+xvec_attr+yvec_attr+"_"+yearstr+".nc")
                    
                    if not os.path.exists(scalar_annual_file): #If it doesn't exist, compute it
                        
                        if scalarECCO:
                            datdirshort, usecompdata = 'Downloads', False
                            
                        elif not scalarECCO:
                            datdirshort, usecompdata = 'computed_monthly', True
                            
                        save_annual_avgs.main(years=[year], field=scalar_attr, datdir=datdirshort, usecompdata=usecompdata, outdir=yearlydatdir)
                    
                    ds_scalar_year = xr.open_mfdataset(scalar_annual_file, engine="scipy")
                    ds_scalar_year.load()

                    #Convert scalar DataSet to useful field
                    lon_centers, lat_centers, lon_edges, lat_edges, scalar_year = ds_to_field(ds_grid, ds_scalar_year.isel(k=k).squeeze(), scalar_attr, latmin, latmax, lonmin, lonmax, resolution)
                    
                    if not os.path.exists(vector_annual_file): #If they don't exist, compute them
                        
                        if vectorECCO:
                            datdirshort, usecompdata = 'Downloads', False
                        
                        elif not vectorECCO:
                            datdirshort, usecompdata = 'computed_monthly', True
                        
                        save_annual_avgs.main(years=[year], field=xvec_attr+yvec_attr, datdir=datdirshort, usecompdata=usecompdata, outdir=yearlydatdir)
                        
                    ds_vector_year = xr.open_mfdataset(vector_annual_file, engine="scipy")
                    ds_vector_year.load()
                    
                    if vectorECCO:
                        vecE, vecN = rotate_vector(ds_grid, ds_vector_year, xvec_attr, yvec_attr)
                        vecE, vecN = vecE.isel(k=k).squeeze(), vecN.isel(k=k).squeeze()
                        
                    elif not vectorECCO:
                        vecE, vecN = rotate_u_g(ds_grid, ds_vector_year[xvec_attr], ds_vector_year[yvec_attr], k)
                    
                    #Plot annual average
                    ArcCir_contourf_quiver(ds_grid, k, [scalar_year], [vecE], [vecN], resolution, cmap, yearstr, lon_centers, lat_centers, scalar_attr, xvec_attr, scalar_bounds=[vmin, vmax], outfile=join(outdir, 'yearly', '{}_k{}_{}.pdf'.format(variables_str, str(k), yearstr)), lats_lons=lats_lons) 
                    
            elif seasonal: #Case where we plot one season per year
                    
                season_months, season_years = get_season_months_and_years(season_start, season_end)
                seas_monthstr = season_months[0] + "-" + season_months[-1] #For titles

                scalar_data_seasons = []
                vecE_data_seasons, vecN_data_seasons = [], []
                
                for i in range(years): #Iterate over specified years
                
                    year = startyr + i
                    yearstr = str(year)
                    endyearstr = str(year + season_years[-1])

                    #Get seasonally-averaged data

                    scalar_seas_file = join(seasonaldatdir, "avg_"+scalar_attr+"_"+season_start+yearstr+"-"+season_end+endyearstr+".nc")
                    vector_seas_file = join(seasonaldatdir, "avg_"+xvec_attr+yvec_attr+"_"+season_start+yearstr+"-"+season_end+endyearstr+".nc")
                        
                    if not os.path.exists(scalar_seas_file): #If it doesn't exist, compute it

                        if scalarECCO:
                            datdirshort, usecompdata = 'Downloads', False

                        elif not scalarECCO:
                            datdirshort, usecompdata = 'computed_monthly', True

                        save_seasonal_avgs.main(field=scalar_attr, years=[year], start_month=season_start, end_month=season_end, usecompdata=usecompdata, datdir=datdirshort, outdir=seasonaldatdir)

                    ds_scalar_seas = xr.open_mfdataset(scalar_seas_file, engine="scipy")
                    ds_scalar_seas.load()
                    
                    #Convert scalar DataSet to useful field
                    lon_centers, lat_centers, lon_edges, lat_edges, scalar_seas = ds_to_field(ds_grid, ds_scalar_seas.isel(k=k), scalar_attr, latmin, latmax, lonmin, lonmax, resolution)
                    
                    #Save seasonal scalar data
                    scalar_data_seasons.append(scalar_seas)
                    
                    if not os.path.exists(vector_seas_file): #If it doesn't exist, compute it

                        if vectorECCO:
                            datdirshort, usecompdata = 'Downloads', False

                        elif not vectorECCO:
                            datdirshort, usecompdata = 'computed_monthly', True

                        save_seasonal_avgs.main(field=xvec_attr+yvec_attr, years=[year], start_month=season_start, end_month=season_end, usecompdata=usecompdata, datdir=datdirshort, outdir=seasonaldatdir)
                            
                    ds_vector_seas = xr.open_mfdataset(vector_seas_file, engine="scipy")
                    ds_vector_seas.load()
                    
                    if vectorECCO:
                        vecE, vecN = rotate_vector(ds_grid, ds_vector_seas, xvec_attr, yvec_attr)
                        vecE, vecN = vecE.isel(k=k).squeeze(), vecN.isel(k=k).squeeze()
                        
                    elif not vectorECCO:
                        vecE, vecN = rotate_u_g(ds_grid, ds_vector_seas[xvec_attr], ds_vector_seas[yvec_attr], k)
                        
                    #Save seasonal vector data
                    
                    vecE_data_seasons.append(vecE)
                    vecN_data_seasons.append(vecN)
                    
                    seas_yearstr = yearstr + "-" + str(year + season_years[-1]) #For titles
                    
                    #Plot seasonal average
                    ArcCir_contourf_quiver(ds_grid, k, [scalar_seas], [vecE], [vecN], resolution, cmap, '{}, {}'.format(seas_monthstr, seas_yearstr), lon_centers, lat_centers, scalar_attr, xvec_attr, scalar_bounds=[vmin, vmax], outfile=join(outdir, 'seasonal', '{}_k{}_{}_{}.pdf'.format(variables_str, str(k), seas_monthstr, seas_yearstr)), lats_lons=lats_lons) 
                    
                #Average over all the seasons and plot
                ArcCir_contourf_quiver(ds_grid, k, scalar_data_seasons, vecE_data_seasons, vecN_data_seasons, resolution, cmap, '{}, {}-{}'.format(seas_monthstr, str(startyr), endyearstr), lon_centers, lat_centers, scalar_attr, xvec_attr, scalar_bounds=[vmin, vmax], outfile=join(outdir, 'seasonal', '{}_k{}_{}_{}-{}.pdf'.format(variables_str, str(k), seas_monthstr, str(startyr), endyearstr)), lats_lons=lats_lons) 
            
##############################

if __name__ == "__main__":
    main()