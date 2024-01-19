"""
Main loop removes saved primary data, if indicated.

Rosalie Cormier, 2024
"""

from functions_remove_data import remove_primary_files

##############################

def main(**kwargs):
    
    if kwargs:
        clear_data_files = kwargs.get('clear_data_files')
    
    if clear_data_files:
        
        if kwargs:
            
            date_strings = kwargs.get('date_strings')
            
            datdir_primary = kwargs.get('datdir_primary')
            
            primary_scalar_fields = kwargs.get('primary_scalar_fields')
            primary_vector_fields = kwargs.get('primary_vector_fields')
            
            secondary_scalar_fields = kwargs.get('secondary_scalar_fields')
            secondary_vector_fields = kwargs.get('secondary_vector_fields')

        #Delete primary scalar data
        for scalar_field_name in primary_scalar_fields:
            remove_primary_files(scalar_field_name, datdir_primary, 
                                 date_strings)

        #Delete primary vector data
        for vector_field_name in primary_vector_fields: 
            remove_primary_files(vector_field_name, datdir_primary, 
                                 date_strings)

        #Identify any secondary fields that required primary data to be saved, 
        #and delete the primary data

        for scalar_field_name in secondary_scalar_fields:
            if scalar_field_name in ['vorticity', 'normal_strain', 
                                     'shear_strain', '2D_div_vel']:
                remove_primary_files('horizontal_vel', datdir_primary, 
                                     date_strings)

        for vector_field_name in secondary_vector_fields:
            if vector_field_name in ['geostrophic_vel', 'Ek_vel']:
                remove_primary_files('density_anom', datdir_primary, 
                                     date_strings)
            if vector_field_name in ['Ek_vel']:
                remove_primary_files('wind_stress', datdir_primary, 
                                     date_strings)

##############################

if __name__ == "__main__":
    main()