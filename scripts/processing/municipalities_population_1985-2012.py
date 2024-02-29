# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 17:30:05 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that formats the data in input file containing the population counts per
municipality for the period 1981 to 2012 into the format:
    year | state | s_code | m_code | m_name | population

which is the format I already have for the period 2013-2020.

UPDATE (21/08/2023): It can also be used with livestock data.

Input files:
    ../../data/raw/BRA_population/municipality_population_1981-2012.csv
    OR
    ../../data/raw/livestock/livestock_data_total.csv
    ../../data/shapes/municipalities_{}.shp
    
Output files: 
    ../../data/processed/BRA_population/municipality_population_{}.csv
    
"""

"""
Imports
"""
import pandas as pd
import geopandas as gpd

import sys
sys.path.append('../utils/')
from municipalities_years_to_shp import municipality_csv_to_shp
corr = municipality_csv_to_shp()

"""
Functions
"""


"""
MAIN
"""

def main():
    
    # --------------------- User inputs ------------------------------------
    
    # File with population counts in the period 1981 - 2012
    data_fp = '../../data/raw/livestock/livestock_data_total.csv' 
    # data_fp = '../../data/raw/BRA_population/municipality_population_1981-2012.csv'
    
    # File to municipalities' shapefiles
    shp_fp = '../../data/shapes/municipalities_shp/municipalities_{}.shp'
    
    # Fiel to save the output
    out_fp = '../../data/processed/livestock/municipalities_livestock_{}.csv'
    # out_fp = '../../data/processed/BRA_population/municipalities_pop_{}.csv'
    
    # Population column name
    col_name = 'animals'
    #col_name = 'population'
    
    start_year = 1985
    
    end_year = 2020
    
    
    # ------------------ Process data --------------------------------------
    
    
    # Read data
    data = pd.read_csv(data_fp)
    
    # Working on one year at a time
    for year in range(start_year, end_year + 1):
        
        # Read the corresponding grid file
        grid_y = gpd.read_file(shp_fp.format(corr[year]))
        
        # Subset relevant columns
        grid_y = grid_y[['state', 's_code','m_code','m_name']]
        
        # Merge the population data
        grid_y = grid_y.merge(data[['m_code', str(year)]], on = 'm_code', how = 'left')
        
        
        # Create column with year
        grid_y['year'] = year
        
        # Renaming population column
        grid_y = grid_y.rename(columns = {str(year): col_name})
        
        # Reorder columns
        grid_y = grid_y[['year', 'state', 's_code', 'm_code', 'm_name', col_name]]
        
        # Save file
        grid_y.to_csv(out_fp.format(year), index = False)
        
        # break
    
if __name__ == '__main__':
    main()

