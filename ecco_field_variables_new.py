"""
Rosalie Cormier, 2023
"""

def get_field_vars(attribute):
    
    monthly_shortnames = {'PHIHYDcR': 'ECCO_L4_DENS_STRAT_PRESS_LLC0090GRID_MONTHLY_V4R4', \
                         'UVEL': 'ECCO_L4_OCEAN_VEL_LLC0090GRID_MONTHLY_V4R4', \
                         'VVEL': 'ECCO_L4_OCEAN_VEL_LLC0090GRID_MONTHLY_V4R4'}
    
    monthly_nc_strings = {'PHIHYDcR': 'OCEAN_DENS_STRAT_PRESS_mon_mean_', \
                         'UVEL': 'OCEAN_VELOCITY_mon_mean_', \
                         'VVEL': 'OCEAN_VELOCITY_mon_mean_'}
    
    return monthly_shortnames[attribute], monthly_nc_strings[attribute]
