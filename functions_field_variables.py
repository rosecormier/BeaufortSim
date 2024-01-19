"""
Contains functions that get each of the following attributes associated with a 
particular field:
    -String ("field variable") used generally in data storage;
    -Strings ("vector comps") associated with x- and y-components of a vector
    variable;
    -Boolean ("field is primary") indicating whether field comes directly from 
    ECCO;
    -Monthly shortname used in monthly-averaged datafile directories;
    -String ("monthly nc string") used in monthly-averaged nc datafiles;
    -Name of colormap to be used in plotting a scalar variable, and whether that
    scalar should be symmetric about zero;
    -Colorbar label to be used in plotting a scalar variable;
    -String to be used in plot titles.
    
Also contains a function that, given lists of scalar and vector fields,
separates each list into a list of primary fields and a list of secondary 
fields.

Rosalie Cormier, 2024
"""

##############################

def get_field_variable(field_name):

    field_variables = {'density_anom': 'RHOAnoma', \
                      'pressure': 'PHIHYDcR', \
                      'vertical_vel': 'WVEL', \
                      'horizontal_vel': 'UVELVVEL', \
                      'vorticity': 'ZETA', \
                      'normal_strain': 'NORMAL', \
                      'shear_strain': 'SHEAR', \
                      '2D_div_vel': 'DIVU', \
                      'geostrophic_vel': 'UGVG', \
                      'Ek_vel': 'UEkVEk', \
                      'wind_stress': 'EXFtauxEXFtauy'}
    
    return field_variables[field_name]

##############################

def get_vector_comps(vector_field_name):

    horizontal_comps = {'horizontal_vel': 'UVEL', \
                       'geostrophic_vel': 'UG', \
                       'Ek_vel': 'UEk', \
                       'wind_stress': 'EXFtaux'}
    vertical_comps = {'horizontal_vel': 'VVEL', \
                     'geostrophic_vel': 'VG', \
                     'Ek_vel': 'VEk', \
                     'wind_stress': 'EXFtauy'}
    
    return (horizontal_comps[vector_field_name], 
            vertical_comps[vector_field_name])

##############################

def field_is_primary(field_name):
    
    field_variable = get_field_variable(field_name)
    
    primary_variables = ['RHOAnoma', 'PHIHYDcR', 'WVEL', 'UVELVVEL', 
                         'EXFtauxEXFtauy']
    secondary_variables = ['ZETA', 'NORMAL', 'SHEAR', 'DIVU', 'UGVG', 'UEkVEk']
    
    if field_variable in primary_variables:
        return True
    elif field_variable in secondary_variables:
        return False

##############################

def get_monthly_shortname(field_variable):
    
    monthly_shortnames = {'RHOAnoma': 
                          'ECCO_L4_DENS_STRAT_PRESS_LLC0090GRID_MONTHLY_V4R4', 
                        'PHIHYDcR': 
                          'ECCO_L4_DENS_STRAT_PRESS_LLC0090GRID_MONTHLY_V4R4',
                        'UVELVVEL': 
                          'ECCO_L4_OCEAN_VEL_LLC0090GRID_MONTHLY_V4R4',
                        'WVEL': 'ECCO_L4_OCEAN_VEL_LLC0090GRID_MONTHLY_V4R4',
                        'UGVG': 'GEOS_VEL_MONTHLY',
                        'ZETA': 'VORTICITY_MONTHLY',
                        'NORMAL': 'STRAIN_MONTHLY',
                        'SHEAR': 'STRAIN_MONTHLY',
                        'EXFtauxEXFtauy': 
                          'ECCO_L4_STRESS_LLC0090GRID_MONTHLY_V4R4',
                        'UEkVEk': 'EK_VEL_MONTHLY',
                        'DIVU': 'DIVU_MONTHLY'}
    
    return monthly_shortnames[field_variable]

##############################
    
def get_monthly_nc_string(field_variable):
    
    monthly_nc_strings = {'RHOAnoma': 'OCEAN_DENS_STRAT_PRESS_mon_mean_',
                          'PHIHYDcR': 'OCEAN_DENS_STRAT_PRESS_mon_mean_',
                         'UVELVVEL': 'OCEAN_VELOCITY_mon_mean_',
                         'WVEL': 'OCEAN_VELOCITY_mon_mean_',
                         'UGVG': 'OCEAN_GEOS_UVEL_mon_mean_',
                         'ZETA': 'OCEAN_VORTICITY_mon_mean_',
                         'NORMAL': 'OCEAN_NORMAL_STRAIN_mon_mean_',
                         'SHEAR': 'OCEAN_SHEAR_STRAIN_mon_mean_',
                         'EXFtauxEXFtauy': 
                          'OCEAN_AND_ICE_SURFACE_STRESS_mon_mean_',
                         'UEkVEk': 'OCEAN_EK_VEL_mon_mean_',
                         'DIVU': 'OCEAN_DIVU_mon_mean_'}
    
    return monthly_nc_strings[field_variable]

##############################

def get_seasonal_shortname(field_variable):
    
    seasonal_shortnames = {'RHOAnoma': 'RHOAnoma_SEASONAL',
                           'PHIHYDcR': 'PHIHYDcR_SEASONAL',
                           'UVELVVEL': 'UVELVVEL_SEASONAL',
                           'WVEL': 'WVEL_SEASONAL', 
                           'UGVG': 'GEOS_VEL_SEASONAL', 
                           'ZETA': 'VORTICITY_SEASONAL', 
                           'NORMAL': 'STRAIN_SEASONAL', 
                           'SHEAR': 'STRAIN_SEASONAL',
                           'EXFtauxEXFtauy': 'EXFtauxEXFtauy_SEASONAL', 
                           'UEkVEk': 'EK_VEL_SEASONAL', 
                           'DIVU': 'DIVU_SEASONAL'}
    
    return seasonal_shortnames[field_variable]

##############################

def get_seasonal_nc_string(field_variable):
    
    seasonal_nc_strings = {'RHOAnoma': 'OCEAN_DENS_STRAT_PRESS_seas_mean_',
                         'PHIHYDcR': 'OCEAN_DENS_STRAT_PRESS_seas_mean_',
                         'UVELVVEL': 'OCEAN_VELOCITY_seas_mean_',
                         'WVEL': 'OCEAN_VELOCITY_seas_mean_',
                         'UGVG': 'OCEAN_GEOS_UVEL_seas_mean_',
                         'ZETA': 'OCEAN_VORTICITY_seas_mean_',
                         'NORMAL': 'OCEAN_NORMAL_STRAIN_seas_mean_',
                         'SHEAR': 'OCEAN_SHEAR_STRAIN_seas_mean_',
                         'EXFtauxEXFtauy': 
                           'OCEAN_AND_ICE_SURFACE_STRESS_seas_mean_',
                         'UEkVEk': 'OCEAN_EK_VEL_seas_mean_',
                         'DIVU': 'OCEAN_DIVU_seas_mean_'}
    
    return seasonal_nc_strings[field_variable]

##############################

def get_cmap_and_symmetry(field_name):
    
    #True iff we want the range of the variable to be centered on zero for 
    #plotting
    symmetry_about_zero = {'density_anom': True,
                          'pressure': False,
                          'vertical_vel': True,
                          'vorticity': True,
                          'normal_strain': True,
                          'shear_strain': True,
                          '2D_div_vel': True}
    
    #MPL colormap name
    cmap = {'density_anom': 'seismic',
            'pressure': 'viridis', 
            'vertical_vel': 'seismic',
            'vorticity': 'seismic',
            'normal_strain': 'seismic',
            'shear_strain': 'seismic', 
            '2D_div_vel': 'seismic'}
    
    return cmap[field_name], symmetry_about_zero[field_name]

##############################

def get_cbar_label(scalar_field_name):

    cbar_label_dict = {'pressure': 
                       r'Hydrostatic pressure anomaly $({m}^2 /{s}^2)$',
                       'density': r'Density anomaly $(kg/{m}^3)$',
                       'vertical_vel': 'Velocity (m/s)',
                       'vorticity': 'Vorticity (1/s)',
                       'normal_strain': r'Normal strain $(1/s^2)$',
                       'shear_strain': r'Shear strain $(1/s^2)$',
                       '2D_div_vel': 'Horizontal velocity divergence (1/s)'}
    
    label = cbar_label_dict[scalar_field_name]
    
    return label

##############################

def get_field_title(field_name):
    
    field_titles = {'RHOAnoma': 'Density Anomaly',
                    'PHIHYDcR': 'Hydrostatic Pressure Anomaly',
                    'WVEL': 'Vertical Velocity',
                    'ZETA': 'Vorticity',
                    'NORMAL': 'Normal Strain',
                    'SHEAR': 'Shear Strain',
                    'DIVU': 'Divergence of Horizontal Velocity',
                    'UVELVVEL': 'Horizontal Velocity',
                    'UGVG': 'Geostrophic Velocity',
                    'EXFtauxEXFtauy': 'Surface Wind-on-Ocean Stress',
                    'UEkVEk': 'Ekman Velocity'}

    field_title = field_titles[get_field_variable(field_name)]
    
    return field_title

##############################

def get_field_lists(scalar_fields, vector_fields):

    primary_scalar_fields, secondary_scalar_fields = [], []

    for scalar_field_name in scalar_fields:
        if field_is_primary(scalar_field_name):
            primary_scalar_fields.append(scalar_field_name)
        elif not field_is_primary(scalar_field_name):
            secondary_scalar_fields.append(scalar_field_name)
        
    primary_vector_fields, secondary_vector_fields = [], []

    if len(vector_fields) != 0:
        for vector_field_name in vector_fields:
            if field_is_primary(vector_field_name):
                primary_vector_fields.append(vector_field_name)
            elif not field_is_primary(vector_field_name):
                secondary_vector_fields.append(vector_field_name)

    return (primary_scalar_fields, secondary_scalar_fields, 
            primary_vector_fields, secondary_vector_fields)