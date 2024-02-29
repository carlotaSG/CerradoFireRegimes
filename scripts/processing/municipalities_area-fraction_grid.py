# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 13:44:33 2022

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script modified from ./_chapter01 on 17/08/2023.

Script calculates, for each cell in a grid, the fraction of area occupied by each
overlapping municipality. It produces a CSV file in the following format:
    year | cell_ID | m_code | frac_area
    
it does so for each shapfile in the input folder (which belong to different years).

Input:
    Grid: ../../data/shapes/cerrado_grid_{}_cropped.shp
    Municipalities: ../../data/shapes/municipalities_shp/municipalities_{}.shp
"""

"""
Imports
"""
import pandas as pd
import geopandas as gpd
from pytictoc import TicToc
t = TicToc()

import sys
sys.path.append('../utils')
import fileList
import calculate_fracAreaClipped as fracArea

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
    grid_id = '30km'
    
    # Grid file
    grid_fp = '../../data/shapes/cerrado_grid_{}_cropped.shp'.format(grid_id)
    
    # # Threshold to select cells
    # area_thresh = 75.0
    
    # include_cells = [567, 568]
    # exclude_cells = [1451]
    
    # Municipalities' shapefiles folder
    in_dir = '../../data/shapes/municipalities_shp/'
    
    # Output folder
    out_dir = '../../data_new_grids/tools/municipalities_areaFraction_grid{}.csv'.format(grid_id)
    
    
    # -------------------------------------------------------------------------
    # Calculations
    
    # Calculating the fraction of area occupied by each municipality in each cell
    
    # We work on municipality SHP individually and consecutively 
    # TIP: if this takes a long time to run, I can run on various nodes for
    # different SHP files
    
    # Read grid file
    grid = gpd.read_file(grid_fp)
    grid = grid.to_crs('epsg:5880')
    
    # # Selecting only those cells whose percentage inside the Cerrado is above a
    # # certain threshold, or the cells that are included, and excluding if necessary
    # grid = grid[(grid['%area_crop'] >= area_thresh) | grid['id'].isin(include_cells)].reset_index(drop = True)
    # grid = grid[~grid['id'].isin(exclude_cells)].reset_index(drop = True)
    
    # print(grid)
    
    # Reading list of municipality SHP files
    m_files = fileList.fileList(in_dir, in_extensions = 'shp')
    t.tic()
    output_data = []

    for file in m_files:
        
        # Read file
        m_data = gpd.read_file(file)
        
        # Getting the year this municipality data belongs to
        year = m_data.loc[0,'year']
        
        print('Working on municipality SHP file for year {}'.format(year))
        
        # Converting CRS to SIRGAS 2000 Brazil Polyconic
        m_data = m_data.to_crs('epsg:5880')
        
        
        # Calling function that calculates the area fractions per cell
        # clipped = fracArea.clip_calculateArea(m_data, grid.loc[0,'geometry'], pol_id = 'm_code', geom_id = grid.loc[0, 'id'])
        
        clipped = fracArea.clip_calculateArea_withPool(m_data, grid, 'm_code','id', num_process = 10)
        # print(clipped)
        
        # Renaming cell_id column (function names it 'geom')
        clipped = clipped.rename(columns = {'geom':'id'})
        # print(clipped)
        
        # Inserting column with the year
        clipped.insert(1, 'year', year)
        
        # print(clipped)
        
        # Append dataframe to list of dataframes
        output_data.append(clipped)
        
        print('Finished working on file {} \n'.format(year))
        # break
    
    # Converting list of dataframes into a single dataframe
    output_data = pd.concat(output_data, ignore_index = True)
    
    
    # Write dataframe to file
    print('Saving file...')
    output_data.to_csv(out_dir, index = False)
    
    
    t.toc('Script ran in')
    
    
    
    
    
    
    
    
    
    
    
    
    
    