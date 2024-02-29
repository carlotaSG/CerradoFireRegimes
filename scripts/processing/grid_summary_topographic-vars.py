# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 11:36:16 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculates the average topographic value per grid cell for each 
topographic variable.
"""

"""
Imports
"""
import geopandas as gpd
from pytictoc import TicToc
t = TicToc()
# import rasterio

import sys
sys.path.append('../utils')
import quantityRasterValue_toPolygon as qty

"""
Functions
"""


"""
MAIN
"""

def main():
    
    # --------------------------- User inputs --------------------------------
    
    # Generic filepath to files with topographic variables
    data_fp = '../../data/processed/SRTM_DEM/SRTM_DEM_merged_{}.tif'
    
    # Topographic variables 
    rough = 'roughness'
    slope = 'slope'
    tri = 'TRI'
    
    list_vars = [rough, slope, tri]
    
    # UPDATE (06/10/2023): including elevation 
    elevation_fp = '../../data/raw/SRTM_DEM/SRTM_DEM_merged.tif'
    
    
    # Grid
    grid_id = '50km'
    grid_fp = '../../data/shapes/cerrado_grid_{}_cropped.shp'.format(grid_id)
    
    # Indications to select cells
    # Threshold to select cells
    area_thresh = 75.0
    # Cells to include/exclude
    include_cells = [567, 568]
    exclude_cells = [1451]
    
    # Output file
    out_fp = '../../data/processed/summary_tables/topography_allVars_grid{}.csv'.format(grid_id)
    
    
    # -------------------- Processing data ----------------------------------


    # TODO: include this after working on datafiles
    # Rading grid 
    grid = gpd.read_file(grid_fp)
    # Selecting only those cells whose percentage inside the Cerrado is above a
    # certain threshold, or the cells that are included, and excluding if necessary
    grid = grid[(grid['%area_crop'] >= area_thresh) | grid['id'].isin(include_cells)].reset_index(drop = True)
    grid = grid[~grid['id'].isin(exclude_cells)].reset_index(drop = True)
    # Subsetting relevant column
    grid = grid[['id']]
    
    # Initializing 
        
    
    # Working on a topographic variable at a time
    for var in list_vars:
        
        t.tic()
        print('-------------------------------------------------------------')
        print('Working on {}...'.format(var))
        
        # Formatting the filepath
        top_fp = data_fp.format(var)
        
        # Calculate the average topographic value per cell in grid
        top_summ = qty.quantityRasterValues_toPolygon(grid_fp, top_fp, nodata = None, quantity = 'mean', num_process = 10)
        
        # Returns a dataframe in the format: id | ... | val
        # Subset relevant columns
        top_summ = top_summ[['id', 'val']]
        
        # Renaming val column with topographic indicator
        top_summ = top_summ.rename(columns = {'val' : var})
        
        # Merging data to grid
        grid = grid.merge(top_summ, on = 'id', how = 'left')
        
        t.toc('{} produced in'.format(var))
        print('')
        
    # UPDATE (06/10/2023): update to include the calculation of average elevation per cell
    t.tic()
    print('-------------------------------------------------------------')
    print('Working on elevation...')
    
    # Formatting the filepath
    top_fp = elevation_fp
    
    # Calculate the average topographic value per cell in grid
    top_summ = qty.quantityRasterValues_toPolygon(grid_fp, top_fp, nodata = -9999, quantity = 'mean', num_process = 10)
    
    # Returns a dataframe in the format: id | ... | val
    # Subset relevant columns
    top_summ = top_summ[['id', 'val']]
    
    # Renaming val column with topographic indicator
    top_summ = top_summ.rename(columns = {'val' : 'elevation'})
    
    # Merging data to grid
    grid = grid.merge(top_summ, on = 'id', how = 'left')
    
    t.toc('Elevation produced in')
    print('')
        
    # Saving output
    grid.to_csv(out_fp, index = False)
        
    
    
if __name__ == '__main__':
    main()
