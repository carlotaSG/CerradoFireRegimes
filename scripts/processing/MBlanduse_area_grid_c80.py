# -*- coding: utf-8 -*-
"""
Created on Wed Dec 22 18:06:50 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that creates a summary CSV file of MapBiomas land use data for a certain
grid. It calculates the following variables for each unit/cell in the grid and
for a specific year. 

It does so using the input file:
    - 'D:/carlota/projects/data/MBlanduse_c06/mapbiomas-brazil-collection-60-cerrado-{}.tif'
    
Variables (-> column_name):
    - anthropic (comprises: Farming (14)) -> 'anthropic'
    - natural (comprises: Forest (1), Wetlands (11), Grassland (12)) -> 'natural'
    - flammable (comprises: Forest (1), Wetlands (11), Grassland (12), 
                 Other non Forest Formations(13), Farming (14)) -> 'flammable'
    - nonflammable (comprises: Salt Flat (32), Rocky Outcrop (29), 
                    Non vegetated area (22), Water(26), Non Observed (27)) -> 'nonflammable'
    - Savanna Formation (4) -> '4'
    - Forest Formation (3) -> '3'
    - Grassland (12) -> '12'
    - Wetlands (11) -> '11'
    - Pasture (15) -> '15'
    - Agriculture (18) -> '18'
    - Forest plantation (9) -> '9'
    - Mosaic of Uses (21) -> '21'
"""

"""
Imports
"""
from pytictoc import TicToc
t = TicToc()

# import geopandas as gpd
import pandas as pd


import sys
sys.path.append('../utils/')

import MB_tools as mbtools
from assignRasterValue_toPolygon import assignRasterValues_toPolygon


"""
Functions
"""
def landuse(data, min_cat, max_cat):
    """
    Function that selects a subset of land-use categories, as well as adding up
    to have the total number of 'natural formation' pixels and 'anthropic formation'
    pixels.

    Parameters
    ----------
    grid : GeoDataFrame
        The GeoDataFrame containing in which we are going to add the percentage 
        area of each land use category in different columns.
    grid_size : integer
        the size in km of the grid's cell, which identifies the grid.
    year : integer
        The year with which we are working.
    min_cat : string
        The label of the land use category with the lowest value.
    max_cat : string
        The label of the land use category with the highest value.

    Returns
    -------
    grid : GeoDataFrame
        The GeoDataFrame containing in which we have added the percentage 
        area of each land use category in different columns.

    """
    
    # Data file
    # data_fp = '../data/processed/grid_{}km/MBfogo_landuse_grid{}km_cerrado_{}.csv'.format(grid_size, grid_size, year)

    # Reading data file
    # data = pd.read_csv(data_fp)
    
    # First, calculate area percentages as the area for a certain category
    # divided by the total area
    data.iloc[:, data.columns.get_loc(min_cat) : data.columns.get_loc(max_cat)] = \
        100 * data.iloc[:, data.columns.get_loc(min_cat) : data.columns.get_loc(max_cat)].div(data['totalArea'], axis = 0)
        
    # Add up total percentage for anthropic categories
    # Update 22/02/2022: Included category 25
    data['anthropic'] = data['14']
    
    # Add up total percentage for natural categories
    data['natural'] = data[['1', '11', '12']].sum(axis = 1)
    
    # Add up total percentage for flammable categories
    data['flammable'] = data[['1', '11', '12', '13','14']].sum(axis = 1)
    
    # Add up total percentage for non-flammable categories
    data['nonflammable'] = data[['32', '29', '22', '26','27']].sum(axis = 1)
    
    # Selecting columns that we are interested in
    data = data[['id', 'anthropic', 'natural', 'flammable', 'nonflammable', '4', '3', '12', '11', '15', '18', '9', '21']]
    
    return data


"""
MAIN
"""

def main():
    
    # User inputs ------------------------------------------------------------
    
    # Choose grid
    grid_size = '50'
    grid_units = 'km'
    grid_id = grid_size + grid_units
    
    # Period over which to work
    # Start year
    start_year = 1985
    # Final year
    end_year = 2022
    
    # Minimum category: the land use category with the smallest number assigned as label
    min_lu = '1'
    # Maximum category: the land use category with the largest number assigned as label
    max_lu = '62'
    # HINT: these two numbers are to be modified depending on the MapBiomas land-use collection being used.
    # MapBiomas Collection 6: min_lu = 1, max_lu = 49
    
    
    # Input files -------------------------------------------------------------
        
    # Land use
    # lu_fp = 'D:/carlota/projects/data/MBlanduse_c06/mapbiomas-brazil-collection-60-cerrado-{}.tif'
    lu_fp = '../../data/raw/MBlulc/mapbiomas-brazil-collection-80-{}.tif'    
    # Grid
    grid_fp = '../../data/shapes/cerrado_grid_{}_cropped.shp'.format(grid_id)
    
    
    # Output files ------------------------------------------------------------
    
    # Land-use in area (km2) per category for all categories
    output_lu_area = '../../data/processed/lulc_data/MBlanduse_grid{}_cerrado_{}.csv'
    
    # Land-use in % area per category for a subset of categories
    output_lu_perc = '../../data/processed/lulc_data/MBlanduse_percArea_subsetLULC_grid{}_{}.csv'
    
    
    # Procesing data ----------------------------------------------------------
    
    for year in range(start_year, end_year + 1):
        
        print('Working on year {}'.format(year))
        
        
        # 1. Count pixels per land use in each grid cell
        # List of numbers within category's range that do not correspond to any 
        # MB land-use category 
        noncat = [0, 2, 7, 8, 16, 17, 28, 34, 37, 38, 42, 43, 44, 45, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61]
    
        t.tic()
        print('Counting pixels per grid cell and category')
        lu = assignRasterValues_toPolygon(grid_fp, lu_fp.format(year), 0, 62, nodata = 99, non_categories = noncat, num_process = 30)
        t.toc('Pixels per land use counted in')
        
        # 2. Discarding columns not of interest
        if grid_size == '50':
            lu = pd.DataFrame(lu.drop(columns = ['withinCerr','%area_crop', 'clipArea', 'geometry','lon','lat']))
        elif grid_size == '30':
            lu = pd.DataFrame(lu.drop(columns = ['withinCerr','clipArea', 'geometry']))
            
        # Filling in MB land use categories that have sub-categories
        lu = mbtools.fill_classes(lu, 8)
        
        
        # 2. Calculating area per cell and land-use category
        for i in range(lu.columns.get_loc('1'), lu.columns.get_loc('62') + 1):
            lu.iloc[:, i] = lu.iloc[:, i] * 30 * 30 * 1e-06
            
        # Saving to file
        lu.to_csv(output_lu_area.format(grid_id, year), index = False)
        
        
        # 3. Percentage area per category for a subset of land-use categories
        
        # Counting number of pixels in each cell as sum of all categories (not falling in noncat)
        # lu['totalArea'] = lu.iloc[:, lu.columns.get_loc('1'):lu.columns.get_loc('49')].sum(axis = 1)
        # Update 22/02/2022: Corrected this to calculate total area as only 
        # the sum of the composite categories
        lu['totalArea'] = lu[['1', '10', '14', '22', '26', '27']].sum(axis = 1)
        


        lu = landuse(lu, min_lu, max_lu)
        
        print('Saving summary file...')

        lu.to_csv(output_lu_perc.format(grid_id, year), index = False)
        t.toc('Year {} summarised in'.format(year))
        
        print('')
        # break
    
if __name__ == '__main__':
    main()
