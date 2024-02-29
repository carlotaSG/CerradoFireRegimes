# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 15:41:33 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculates the area in km2 per cell in grid occupied by a protected
area for each protected area file.


"""

"""
Imports
"""
import pandas as pd
import geopandas as gpd
from pytictoc import TicToc
t = TicToc()

import sys
sys.path.append('../utils/')
import calculate_fracAreaClipped as fracArea

"""
Functions
"""



"""
MAIN
"""
def main():
    
    # -------------------- User inputs ------------------------------------
    
    # Files to work with
    in_dir = '../../data/shapes/protected_areas/'
    
    ti_fp = 'pa_terras-indigenas.shp'
    ti_pi_fp = 'pa_terras-indigenas_protecao-integral.shp'
    us_fp = 'pa_uso-sustentavel.shp'
    all_fp = 'pa_all-areas.shp'

    # List of files and the words to use to store their areas in file
    list_files = list(zip(
        [ti_fp, ti_pi_fp, us_fp, all_fp],
        ['ti', 'ti_pi', 'us', 'all']
        ))
    
    # The grid
    grid_id = '50km'
    
    # Threshold to select cells
    area_thresh = 75.0
    
    include_cells = [567, 568]
    exclude_cells = [1451]
    
    grid_fp = '../../data/shapes/cerrado_grid_{}_cropped.shp'.format(grid_id)
    
    # Output filepath
    output_fp = '../../data/processed/summary_tables/protected_areas_{}.csv'.format(grid_id)
    
    
    
    # -------------------- Processing -------------------------------------
    
    # Reading the grid and subsetting relevant cells
    grid = gpd.read_file(grid_fp)
    grid = grid[(grid['%area_crop'] >= area_thresh) | grid['id'].isin(include_cells)]
    grid = grid[~grid['id'].isin(exclude_cells)].reset_index(drop = True)
    
    # Adjusting CRS to metered one
    grid = grid.to_crs('epsg:5880')
    
    # Subsetting relevant columns only
    grid = grid[['id', 'geometry']]
    
    # Working on a protected areas file at a time
    for (file_fp, file_id) in list_files:
        
        # Read protected areas' polygons
        pa_data = gpd.read_file(in_dir + file_fp)
        
        # Converting projection to metered one
        pa_data = pa_data.to_crs('epsg:5880')
        
        # Calculating the area fraction of each protected area inside each grid cell
        pa_areas = fracArea.clip_calculateArea_withPool(pa_data, grid, 'pa_id','id', num_process = 10)
        
        # Reanming columns
        pa_areas = pa_areas.rename(columns = {'geom': 'id', 'frac_area':file_id + '_area'})
        
        # Summarise to have the total protected area fraction per cell
        pa_areas = pa_areas[['id', file_id + '_area']].groupby(by = 'id').sum().reset_index()
        
        # Adding information to grid
        grid = grid.merge(pa_areas, on = 'id', how = 'inner')
        
    # Transform to plain data frame
    grid = grid.drop(columns = 'geometry')
    
    # Save file 
    grid.to_csv(output_fp, index = False)
    
    
    

if __name__ == '__main__':
    main()

