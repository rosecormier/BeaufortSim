"""
Contains functions that get attributes associated with a particular field.
    
Also contains a function that, given lists of scalar and vector fields,
separates each list into a list of primary fields and a list of secondary 
fields.

R. Cormier, 2024
"""

##############################

def get_variable_str(field_name):
    field_variables = {'density_anom': 'RHOAnoma', \
                      'pressure': 'PHIHYDcR', \
                      'vertical_vel': 'WVEL', \
                      'horizontal_vel': 'UVELVVEL', \
                      'uvelwvel': 'UVELWVEL', \
                      'vorticity': 'ZETA', \
                      'normal_strain': 'NORMAL', \
                      'shear_strain': 'SHEAR', \
                      '2D_div_vel': 'DIVU', \
                      'geostrophic_vel': 'UGVG', \
                      'Ek_vel': 'UEkVEk', \
                      'wind_stress': 'EXFtauxEXFtauy', \
                      'ssh': 'SSH', \
                      'potential_temp': 'THETA', \
                      'ice_thickness': 'SIheff', \
                      'T_zflux': 'T_ZFLUX'}
    return field_variables[field_name]

##############################

def get_vector_comps(vector_field_name):

    horizontal_comps = {'horizontal_vel': 'UVEL', \
                       'uvelwvel': 'UVEL', \
                       'geostrophic_vel': 'UG', \
                       'Ek_vel': 'UEk', \
                       'wind_stress': 'EXFtaux'}
    vertical_comps = {'horizontal_vel': 'VVEL', \
                     'uvelwvel': 'WVEL', \
                     'geostrophic_vel': 'VG', \
                     'Ek_vel': 'VEk', \
                     'wind_stress': 'EXFtauy'}
    return (horizontal_comps[vector_field_name], 
                vertical_comps[vector_field_name])

##############################

def field_is_primary(field_name):
    
    variable_str = get_variable_str(field_name)
    
    primary_variables = ['RHOAnoma', 
                         'PHIHYDcR', 
                         'WVEL', 
                         'UVELVVEL', 
                         'UVELWVEL',
                         'EXFtauxEXFtauy', 
                         'SSH', 
                         'THETA', 
                         'SIheff']
    secondary_variables = ['ZETA', 
                           'NORMAL', 
                           'SHEAR', 
                           'DIVU', 
                           'UGVG', 
                           'UEkVEk', 
                           'T_ZFLUX']
    if variable_str in primary_variables:
        return True
    elif variable_str in secondary_variables:
        return False

##############################

def get_monthly_shortname(field_variable):
    monthly_shortnames = {'RHOAnoma': 
             'ECCO_L4_DENS_STRAT_PRESS_LLC0090GRID_MONTHLY_V4R4', 
        'PHIHYDcR': 'ECCO_L4_DENS_STRAT_PRESS_LLC0090GRID_MONTHLY_V4R4',
        'UVELVVEL': 'ECCO_L4_OCEAN_VEL_LLC0090GRID_MONTHLY_V4R4',
        'UVELWVEL': 'ECCO_L4_OCEAN_VEL_LLC0090GRID_MONTHLY_V4R4',
        'WVEL': 'ECCO_L4_OCEAN_VEL_LLC0090GRID_MONTHLY_V4R4',
        'UGVG': 'GEOS_VEL_MONTHLY',
        'ZETA': 'VORTICITY_MONTHLY',
        'NORMAL': 'STRAIN_MONTHLY',
        'SHEAR': 'STRAIN_MONTHLY',
        'EXFtauxEXFtauy': 'ECCO_L4_STRESS_LLC0090GRID_MONTHLY_V4R4',
        'UEkVEk': 'EK_VEL_MONTHLY',
        'DIVU': 'DIVU_MONTHLY', 
        'SSH': 'ECCO_L4_SSH_LLC0090GRID_MONTHLY_V4R4',
        'THETA': 'ECCO_L4_TEMP_SALINITY_LLC0090GRID_MONTHLY_V4R4',
        'SIheff': 'ECCO_L4_SEA_ICE_CONC_THICKNESS_LLC0090GRID_MONTHLY_V4R4',
        'T_ZFLUX': 'T_ZFLUX_MONTHLY'}
    return monthly_shortnames[field_variable]

##############################
    
def get_monthly_nc_string(field_variable):
    monthly_nc_strings = {'RHOAnoma': 'OCEAN_DENS_STRAT_PRESS_mon_mean_',
        'PHIHYDcR': 'OCEAN_DENS_STRAT_PRESS_mon_mean_',
        'UVELVVEL': 'OCEAN_VELOCITY_mon_mean_',
        'UVELWVEL': 'OCEAN_VELOCITY_mon_mean_',
        'WVEL': 'OCEAN_VELOCITY_mon_mean_',
        'UGVG': 'OCEAN_GEOS_UVEL_mon_mean_',
        'ZETA': 'OCEAN_VORTICITY_mon_mean_',
        'NORMAL': 'OCEAN_NORMAL_STRAIN_mon_mean_',
        'SHEAR': 'OCEAN_SHEAR_STRAIN_mon_mean_',
        'EXFtauxEXFtauy': 'OCEAN_AND_ICE_SURFACE_STRESS_mon_mean_',
        'UEkVEk': 'OCEAN_EK_VEL_mon_mean_',
        'DIVU': 'OCEAN_DIVU_mon_mean_',
        'SSH': 'SEA_SURFACE_HEIGHT_mon_mean_',
        'THETA': 'OCEAN_TEMPERATURE_SALINITY_mon_mean_',
        'SIheff': 'SEA_ICE_CONC_THICKNESS_mon_mean_',
        'T_ZFLUX': 'OCEAN_T_ZFLUX_mon_mean_'}
    return monthly_nc_strings[field_variable]

##############################

def get_field_variable(field_variable):
    return (get_monthly_shortname[field_variable], 
                get_monthly_nc_strings[field_variable])
    
##############################

def get_seasonal_shortname(field_variable):
    
    seasonal_shortnames = {'RHOAnoma': 'RHOAnoma_SEASONAL',
                           'PHIHYDcR': 'PHIHYDcR_SEASONAL',
                           'UVELVVEL': 'UVELVVEL_SEASONAL',
                           'UVELWVEL': 'UVELWVEL_SEASONAL',
                           'WVEL': 'WVEL_SEASONAL', 
                           'UGVG': 'GEOS_VEL_SEASONAL', 
                           'ZETA': 'VORTICITY_SEASONAL', 
                           'NORMAL': 'STRAIN_SEASONAL', 
                           'SHEAR': 'STRAIN_SEASONAL',
                           'EXFtauxEXFtauy': 'EXFtauxEXFtauy_SEASONAL', 
                           'UEkVEk': 'EK_VEL_SEASONAL', 
                           'DIVU': 'DIVU_SEASONAL',
                           'SSH': 'SSH_SEASONAL',
                           'THETA': 'THETA_SEASONAL',
                           'SIheff': 'SIheff_SEASONAL',
                           'T_ZFLUX': 'T_ZFLUX_SEASONAL'}
    
    return seasonal_shortnames[field_variable]

##############################

def get_seasonal_nc_string(field_variable):
    seasonal_nc_strings = {'RHOAnoma': 'OCEAN_DENS_STRAT_PRESS_seas_mean_',
        'PHIHYDcR': 'OCEAN_DENS_STRAT_PRESS_seas_mean_',
        'UVELVVEL': 'OCEAN_VELOCITY_seas_mean_',
        'UVELWVEL': 'OCEAN_VELOCITY_seas_mean_',
        'WVEL': 'OCEAN_VELOCITY_seas_mean_',
        'UGVG': 'OCEAN_GEOS_UVEL_seas_mean_',
        'ZETA': 'OCEAN_VORTICITY_seas_mean_',
        'NORMAL': 'OCEAN_NORMAL_STRAIN_seas_mean_',
        'SHEAR': 'OCEAN_SHEAR_STRAIN_seas_mean_',
        'EXFtauxEXFtauy': 'OCEAN_AND_ICE_SURFACE_STRESS_seas_mean_',
        'UEkVEk': 'OCEAN_EK_VEL_seas_mean_',
        'DIVU': 'OCEAN_DIVU_seas_mean_',
        'SSH': 'SSH_seas_mean_',
        'THETA': 'THETA_seas_mean_',
        'SIheff': 'SIheff_seas_mean_',
        'T_ZFLUX': 'T_ZFLUX_seas_mean_'}
    return seasonal_nc_strings[field_variable]

##############################

def get_cmap_and_symmetry(field_name):
    #True iff we want the range to be centered on zero for plotting
    symmetry_about_zero = {'density_anom': False,
                          'pressure': False,
                          'vertical_vel': True,
                          'vorticity': True,
                          'normal_strain': True,
                          'shear_strain': True,
                          '2D_div_vel': True,
                          'ssh': False,
                          'ice_thickness': False,
                          'T_zflux': True}
    #MPL colormap name
    cmap = {'density_anom': 'Spectral_r',
            'pressure': 'viridis', 
            'vertical_vel': 'seismic',
            'vorticity': 'seismic',
            'normal_strain': 'seismic',
            'shear_strain': 'seismic', 
            '2D_div_vel': 'seismic',
            'ssh': 'plasma',
            'ice_thickness': 'Blues_r',
            'T_zflux': 'coolwarm'}
    return cmap[field_name], symmetry_about_zero[field_name]

##############################

def get_cbar_label(scalar_field_name):
    cbar_label_dict = {'pressure': 
            r'Hydrostatic pressure anomaly $({m}^2 /{s}^2)$',
        'density_anom': r'Density anomaly $(kg/{m}^3)$',
        'vertical_vel': 'Velocity (m/s)',
        'vorticity': 'Vorticity (1/s)',
        'normal_strain': r'Normal strain $(1/s^2)$',
        'shear_strain': r'Shear strain $(1/s^2)$',
        '2D_div_vel': 'Horizontal velocity divergence (1/s)',
        'ssh': 'Height (m)',
        'ice_thickness': 'Thickness (m)',
        'T_zflux': r'Flux (m${}^2$K/s)'}
    return cbar_label_dict[scalar_field_name]

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
                    'UEkVEk': 'Ekman Velocity',
                    'SSH': 'Dynamic Sea-Surface Height Anomaly',
                    'SIheff': 'Sea-Ice Thickness',
                    'T_ZFLUX': 'Vertical Potential-Temperature Flux'}
    field_title = field_titles[get_variable_str(field_name)]
    
    #Whether the variable is (can be) defined at multiple depth indices
    multiple_depths_dict = {'RHOAnoma': True,
                            'PHIHYDcR': True,
                            'WVEL': True,
                            'ZETA': True,
                            'NORMAL': True,
                            'SHEAR': True,
                            'DIVU': True,
                            'UVELVVEL': True,
                            'UGVG': True,
                            'EXFtauxEXFtauy': False,
                            'UEkVEk': True,
                            'SSH': False,
                            'SIheff': False,
                            'T_ZFLUX': False}
    multiple_depths = multiple_depths_dict[get_variable_str(field_name)]
    
    return field_title, multiple_depths

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