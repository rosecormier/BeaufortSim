"""
Vorticity-related functions, including Ro_l and OW visualization.

Rosalie Cormier, 2023
"""

import xarray as xr
import ecco_v4_py as ecco

from os.path import join

from ecco_field_variables import get_field_vars
from ecco_general import ecco_resample, load_dataset
from geostrophic_functions import comp_f
from ecco_visualization import ArcCir_pcolormesh

import load_ECCO_dataset

##############################

def comp_vorticity(grid_llc, u_mean, v_mean, dx, dy, cell_area):
    return (-grid_llc.diff(u_mean * dx, 'Y') + grid_llc.diff(v_mean * dy, 'X')).squeeze() / cell_area

##############################

def comp_normal_strain(grid_llc, u_mean, v_mean, dx, dy, cell_area):
    return (grid_llc.diff(u_mean * dy, 'X') - grid_llc.diff(v_mean * dx, 'Y')) / cell_area

##############################

def comp_shear_strain(grid_llc, u_mean, v_mean, dx, dy, cell_area):
    return (grid_llc.diff(v_mean * dy, 'X') + grid_llc.diff(u_mean * dx, 'Y')) / cell_area

##############################

def comp_total_strain(ds_grid, normal_strain, shear_strain, latmin, latmax, lonmin, lonmax, resolution):

    normal_strain_field = ecco_resample(ds_grid, normal_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
    shear_strain_field = ecco_resample(ds_grid, shear_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
    
    return normal_strain_field + shear_strain_field

##############################

def comp_OkuboWeiss(omega, normal_strain, shear_strain):
    return normal_strain**2 + shear_strain**2 - omega**2

##############################

def comp_local_Ro(omega, y):
    
    f = comp_f(y)
    
    return abs(omega) / f

##############################

def plot_monthly_Ro_l(Ro_l_list, zeta_field, lon_centers, lat_centers, seasonal, outdir, k, monthstr, yearstr, ds_grid, resolution, datestr, lats_lons, scalar_bounds=[1e-4, 1e-2]):

    """
    Computes and plots local Rossby number corresponding to monthly vorticity field.
    """
    
    Ro_l = comp_local_Ro(zeta_field, lat_centers) #Compute Ro_l
            
    #Define output file name
        
    if not seasonal:
        Ro_l_outfile = join(outdir, 'monthly', 'localRo_k{}_{}{}.pdf'.format(str(k), monthstr, yearstr))
    elif seasonal:
        Ro_l_outfile = join(outdir, 'seasonal', 'localRo_k{}_{}_{}.pdf'.format(str(k), seas_monthstr, seas_yearstr))

    #Plot Ro_l
    ArcCir_pcolormesh(ds_grid, [Ro_l], resolution, 'Reds', lon_centers, lat_centers, None, datestr, 'Ro_l', scalar_bounds=scalar_bounds, k_plot=k, extend='both', logscale=True, outfile=Ro_l_outfile, lats_lons=lats_lons)

    #Save Ro_l data and return it
    
    Ro_l_list.append(Ro_l)
    return Ro_l_list
    
##############################

def plot_monthly_OW(OW_list, zeta_field, seasonal, yearstr, year, k, datdirname, ds_grid, lon_centers, lat_centers, latmin, latmax, lonmin, lonmax, resolution, datestr, lats_lons, monthstr=None, datdir=None, season_start=None, season_end=None, endyearstr=None, seas_monthstr=None, seas_yearstr=None, seasonaldatdir=None, scalar_bounds=[-1e-14, 1e-14]):
    
    """
    Computes and plots Okubo-Weiss parameter corresponding to monthly velocity and vorticity profiles.
    """
    
    if not seasonal:
        
        #Get monthly velocity data
        
        vel_monthly_shortname, vel_monthly_nc_str = get_field_vars('UVELVVEL')
        vel_file = join(datdir, vel_monthly_shortname, vel_monthly_nc_str+yearstr+"-"+monthstr+"_ECCO_V4r4_native_llc0090.nc")
        ds_vel = load_dataset(vel_file)
        (ds_vel['UVEL']).data, (ds_vel['VVEL']).data = (ds_vel['UVEL']).values, (ds_vel['VVEL']).values
        
        #Define output file name
        OW_outfile = join(outdir, 'monthly', 'OW_k{}_{}{}.pdf'.format(str(k), monthstr, yearstr))
        
    elif seasonal:
        
        #Get seasonal velocity data
        
        vel_seas_file = join(seasonaldatdir, "avg_UVELVVEL_"+season_start+yearstr+"-"+season_end+endyearstr+".nc")
        ds_vel = xr.open_mfdataset(vel_seas_file, engine="scipy")
        ds_vel.load()
        
        #Define output file name
        OW_outfile = join(outdir, 'seasonal', 'OW_k{}_{}_{}.pdf'.format(str(k), seas_monthstr, seas_yearstr))
        
    #Compute strain terms
    
    xgcm_grid = ecco.get_llc_grid(ds_grid)
    normal_strain = comp_normal_strain(xgcm_grid, ds_vel['UVEL'], ds_vel['VVEL'], ds_grid.dxG, ds_grid.dyG, ds_grid.rA).isel(k=k).squeeze()
    shear_strain = comp_shear_strain(xgcm_grid, ds_vel['UVEL'], ds_vel['VVEL'], ds_grid.dxC, ds_grid.dyC, ds_grid.rAz).isel(k=k).squeeze()
    normal_strain= ecco_resample(ds_grid, normal_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
    shear_strain = ecco_resample(ds_grid, shear_strain, latmin, latmax, lonmin, lonmax, resolution)[4]
    
    OW = comp_OkuboWeiss(zeta_field, normal_strain, shear_strain) #Compute OW 
    
    #Plot OW
    ArcCir_pcolormesh(ds_grid, [OW], resolution, 'seismic', lon_centers, lat_centers, None, datestr, 'OW', scalar_bounds=scalar_bounds, k_plot=k, extend='both', outfile=OW_outfile, lats_lons=lats_lons)

    OW_list.append(OW)
    return OW_list