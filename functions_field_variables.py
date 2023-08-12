"""
Rosalie Cormier, 2023
"""

def get_field_vars(attribute):
    
    monthly_shortnames = {'PHIHYDcR': 'ECCO_L4_DENS_STRAT_PRESS_LLC0090GRID_MONTHLY_V4R4', \
                         'UVELVVEL': 'ECCO_L4_OCEAN_VEL_LLC0090GRID_MONTHLY_V4R4', \
                         'WVEL': 'ECCO_L4_OCEAN_VEL_LLC0090GRID_MONTHLY_V4R4', \
                         'UGVG': 'GEOS_VEL_MONTHLY', \
                         'ZETA': 'VORTICITY_MONTHLY', \
                         'NORMAL': 'STRAIN_MONTHLY', \
                         'SHEAR': 'STRAIN_MONTHLY', \
                         'EXFtauxEXFtauy': 'ECCO_L4_STRESS_LLC0090GRID_MONTHLY_V4R4', \
                         'UEkVEk': 'EK_VEL_MONTHLY', \
                         'DIVU': 'DIVU_MONTHLY', \
                         'DIVUEk': 'DIVUEk_MONTHLY', \
                         'SIheff': 'ECCO_L4_SEA_ICE_CONC_THICKNESS_LLC0090GRID_MONTHLY_V4R4', \
                         'EXFuwindEXFvwind': 'ECCO_L4_ATM_STATE_LLC0090GRID_MONTHLY_V4R4', \
                         'SALT': 'ECCO_L4_TEMP_SALINITY_LLC0090GRID_MONTHLY_V4R4'}
    
    monthly_nc_strings = {'PHIHYDcR': 'OCEAN_DENS_STRAT_PRESS_mon_mean_', \
                         'UVELVVEL': 'OCEAN_VELOCITY_mon_mean_', \
                         'WVEL': 'OCEAN_VELOCITY_mon_mean_', \
                         'UGVG': 'OCEAN_GEOS_UVEL_mon_mean_', \
                         'ZETA': 'OCEAN_VORTICITY_mon_mean_', \
                         'NORMAL': 'OCEAN_NORMAL_STRAIN_mon_mean_', \
                         'SHEAR': 'OCEAN_SHEAR_STRAIN_mon_mean_', \
                         'EXFtauxEXFtauy': 'OCEAN_AND_ICE_SURFACE_STRESS_mon_mean_', \
                         'UEkVEk': 'OCEAN_EK_VEL_mon_mean_', \
                         'DIVU': 'OCEAN_DIVU_mon_mean_', \
                         'DIVUEk': 'OCEAN_DIVUEk_mon_mean_', \
                         'SIheff': 'SEA_ICE_CONC_THICKNESS_mon_mean_', \
                         'EXFuwindEXFvwind': 'ATM_SURFACE_TEMP_HUM_WIND_PRES_mon_mean_', \
                         'SALT': 'OCEAN_TEMPERATURE_SALINITY_mon_mean_'}
    
    return monthly_shortnames[attribute], monthly_nc_strings[attribute]

def get_variable_str(attribute, geostrophic=False):
    
    variables = {'PHIHYDcR': 'p_anom', \
                'UVELVVEL': 'u', \
                'WVEL': 'w', \
                'ZETA': 'zeta', \
                'UGVG': 'u_g', \
                'EXFtauxEXFtauy': 'tau', \
                'UEkVEk': 'u_Ek', \
                'DIVU': 'div_u',
                'DIVUEk': 'div_u_Ek', \
                'SIheff': 'icethickness'}

    variable_string = variables[attribute]
        
    return variable_string