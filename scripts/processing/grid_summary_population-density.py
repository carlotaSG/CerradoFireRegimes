# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 16:16:45 2022

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

UPDATE 18/08/2023 to read the correspondence from a script.
UPDATE 22/08/2023 to calculate livestock density - data stored in same format as 
                  human population density - per grid cell and year.

Script that calcualtes the area-weighted average population density per cell in
grid and year.

It reads the population density at municipality level from year files, and then
reads the cells with the corresponding municipalities overlapping it each year along with
each municiaplity's fraction of cell's area.

With these two pieces of information, it calculates the average population density
per cell and year.

It takes into account the fact that municiaplities' areas were only available for a handful
of years, which may correspond to a year in the municiaplity population files, 
or a few years.

Input files:
    Municipalities' area fraction per cell: '../../data/intermediate/municipalities_areaFraction_grid{grid_id}.csv'
    Municipalities' population density: '../../data/processed/BRA_population/municipalities_pop_{year}.csv'
    Municipalities' livestock density: '../../data/processed/livestock/municipalities_livestock_{year}.csv'
    
Output file:
    ../../data/processed/summary_tables/population_density_grid{grid_id}_{year}.csv
    ../../data/processed/summary_tables/livestock_density_grid{grid_id}_{year}.csv
"""


"""
Imports
"""
import pandas as pd

import sys
sys.path.append('../../scripts')
# import fileList
# import calculate_fracAreaClipped as fracArea

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
    
    # Grid selection
    grid_id = '50km'
    
    # User selects the years 
    start_year, end_year = 1985, 2020
    
    # Loading the correspondence between the municipality fraction area per cell
    # and the population density per municipality
    correspondence = municipality_csv_to_shp()
    
    # Input files
    # pop_dens_fp = '../../data/processed/BRA_population/municipalities_pop_{}.csv'
    pop_dens_fp = '../../data/processed/livestock/municipalities_livestock_{}.csv'
    frac_area_fp = '../../data/tools/municipalities_areaFraction_grid{}.csv'.format(grid_id)
    
    # Output file
    # out_fp = '../../data/processed/summary_tables/population_density_grid{}_{}.csv'
    out_fp = '../../data/processed/summary_tables/livestock_density_grid{}_{}.csv'
    
    
    # -------------------------------------------------------------------------
    # Area-weighted population density
    
    # Read fraction of area of municiaplities inside wach cell for all years
    frac_area = pd.read_csv(frac_area_fp)
    
    # Column to read municipality density from
    dens_col = 'anim_density' # 'pop_density'
    
    # Column to store the area weighted density
    areaw_dens = 'areaw_dens' # 'areaw_popdensity'
    
    
    for year in range(start_year, end_year + 1):
        
        print('Working on year {}'.format(year))
        
        # Select the "fraction area" year that corresponds to the time series
        # year
        frac_area_year = correspondence[year]
        
        # print(frac_area_year)
        
        # Read the population per municipality for this year
        pop_dens = pd.read_csv(pop_dens_fp.format(year))

        # Using the "area fraction" year to subset the fraction areas per cell
        # we are interested in
        frac_area_year = frac_area[frac_area['year'] == frac_area_year]
        
        # print(pop_dens)
        # print(frac_area_year)
        
        # Substituting the "area fraction" year by the time series year
        frac_area_year['year'] = year
        
        # merge to it the population density for each municipality
        # Merging the population density for each municipality
        frac_area_year = frac_area_year.merge(pop_dens[['m_code', dens_col]], on = ['m_code'], how = 'left')
        
        # Creating column that contains the product between the area weights and the population density
        frac_area_year[areaw_dens] = frac_area_year['frac_area'] * frac_area_year[dens_col]


        # Calculate the area-weighted population density
        frac_area_year = frac_area_year[['id','year', areaw_dens]].groupby(by = ['id','year']).sum().reset_index()
        
        
        # Now we have the dataframe that we want.
        # Save file for this year
        frac_area_year.to_csv(out_fp.format(grid_id, year), index = False)
        
        # Freeing memory
        del frac_area_year
        # break