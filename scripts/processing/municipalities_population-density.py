# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 11:52:36 2022

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script modified from ./_chapter01 on 17/08/2023 to caluclate the municipalities'
population densities for the years 1992 to 1995.

Script that calculates the population density per year in a certain period
per municipality. The output files have the following information:
    year | state | s_code | m_code | m_name | population | area_km2 | pop_density

It does so by reading the areas from the municipalities' shapefiles for the different
years. Note there are no shapefiles for each year, the correspondence is found in 
the function municipality_csv_to_shp.
    
Input files:
    Population counts: ../../data/processed/BRA_population/municipalities_pop_<year>.csv
    OR
    Animal counts: ../../data/processed/livestock/municipalities_livestock_<year>.csv
    Municipality geometry: ../../data/shapes/municipalities_shp/municipalities_<year>.shp
    
Output files: (input is overwritten)
    ../data/processed/BRA_population/municiaplities_pop_<year>.csv
    
"""

"""
Imports
"""
import pandas as pd
import geopandas as gpd

import sys
sys.path.append('../utils/')
from municipalities_years_to_shp import municipality_csv_to_shp
csv_to_shp = municipality_csv_to_shp()

"""
Functions
"""


"""
MAIN
"""
if __name__ == '__main__':
    
    
    # -------------------------------------------------------------------------
    # User inputs
    
    # The user must declare what years it want to work on
    # Start and end years
    start_year, end_year = 1985, 2020
    
    # Declare input files
    # csv_fp = '../../data/processed/BRA_population/municipalities_pop_{}.csv'
    csv_fp = '../../data/processed/livestock/municipalities_livestock_{}.csv'
    shp_fp = '../../data/shapes/municipalities_shp/municipalities_{}.shp'
    
    # Output files
    out_fp = csv_fp
    
    # Column name for which to calculate the density
    col_name = 'animals' # 'population'
    dens_name = 'anim_density' # 'pop_density'
    
    # -------------------------------------------------------------------------
    # Calculating population density per municipality and year
    
    # Loading correspondance between CSV and SHP years
    csv_to_shp = municipality_csv_to_shp()
    
    # We operate in each year separately
    for year in range(start_year, end_year + 1):
        
        print(year)
        
        # Get the corresponding SHP year
        shp_year = csv_to_shp[year]
        
        # Reading CSV and SHP files
        csv_data = pd.read_csv(csv_fp.format(year))
        shp_data = gpd.read_file(shp_fp.format(shp_year))
        
        # Selecting only relevant information from SHP
        shp_data = shp_data[['m_code','area_km2']]
        
        # Merging information to CSV
        csv_data = csv_data.merge(shp_data, on = ['m_code'], how = 'outer')
        
        # Freeing memory
        del shp_data
        
        # print(csv_data.dtypes)
        
        # Calculate population density
        csv_data[dens_name] = csv_data[col_name] / csv_data['area_km2']
        
        
        # Finally, save data back
        csv_data.to_csv(csv_fp.format(year), index = False)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    