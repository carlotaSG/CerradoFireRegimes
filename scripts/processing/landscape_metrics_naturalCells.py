# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 10:43:10 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculates a series of class-level metrics for the land-use types
considered as natural, anthropic flammable, and non-flammable as a means to 
quantify landscape's spatial heterogeneity. 

Re-categorization:
    natural_cats = [3,4,5,49,11,12,32,29,13]
    anthropic_flammable_cats = [15, 19,39,20,40,41,36,46,47,48,9,21]
    nonflammable_cats = [23,24,30,25,33,31]

It does so using the package "pylandstats", which implements in Python the 
metrics described in FRAGSTATS.

It calculates the class metrics for each cell and year in the period.


"""

"""
Imports
"""
import pandas as pd
import numpy as np
import geopandas as gpd
import rasterio
from rasterio.mask import mask
import pylandstats as pls
from multiprocessing import Pool

from pytictoc import TicToc
t = TicToc()

"""
Functions
"""
def read_grid_rasterCRS(grid_fp, raster_fp, year, thresh = 75.0):
    
    # Get the raster's CRS
    with rasterio.open(raster_fp.format(year)) as src:
        raster_crs = src.crs
        
    
    # Read the grid
    grid = gpd.read_file(grid_fp)
    
    # Transform to rsater's CRS
    grid = grid.to_crs(raster_crs)
    
    # Filter only cell with x% of the area in Cerrado
    grid = grid[grid['%area_crop'] >= thresh]
    
    return(grid)


def calculate_classMetrics(cell_id, polygon, raster_fp, year, list_metrics, list_metrics_kws, no_data):
    

    with rasterio.open(raster_fp.format(year)) as raster:
        
        cell_array, cell_transform = mask(dataset = raster, shapes=[polygon], crop=True, nodata = no_data)
        
        # Reducing the number of dimensions (as there is an empty third dimension)
        cell_array = np.squeeze(cell_array)
        
    natural_cats = [3,4,5,49,11,12,32,29,13]
    anthropic_flammable_cats = [15, 19,39,20,40,41,36,46,47,48,9,21]
    nonflammable_cats = [23,24,30,25,33,31]
        
    # Relabelling together the natural vegetation categories
    cell_array = np.where(np.isin(cell_array, natural_cats), 61, cell_array)
    # Relabelling together the flammable non-natural categories
    cell_array = np.where(np.isin(cell_array, anthropic_flammable_cats), 62, cell_array)
    
    # Relabelling as no data the remaining categories
    cell_array = np.where(np.isin(cell_array, nonflammable_cats), 63, cell_array)
    
    # cell_array.to_csv('../test_{}_{}.csv'.format(year, cell_id), index = False)
        

    # Creating the Landscape object
    landsc = pls.Landscape(cell_array, res = (30, 30), nodata = no_data, 
                           transform = cell_transform)

    # Calculate the metrics
    metrics_df = landsc.compute_class_metrics_df(metrics = list_metrics)
    
    # Casting inde x(the class identifier) as column
    metrics_df = metrics_df.reset_index()
    
    # Relabelling the class categories
    metrics_df['class_val'] = metrics_df['class_val'].map({61: 'natural', 62: 'anthropic', 63: 'non-flammable'})
    
    # Including column with year and cell_id
    metrics_df.insert(loc = 0, column = 'year', value = year)
    metrics_df.insert(loc = 0, column = 'id', value = cell_id)
    
    return metrics_df



def calculate_classMetrics_withPool(grid, raster_fp, year, list_metrics, list_metrics_kws, no_data, num_process = 10):
    
    # Preparing input for Pool
    n_cells = len(grid)
    poolInput = list(zip(grid['id'].tolist(),
                         grid['geometry'].tolist(),
                         [raster_fp]*n_cells,
                         [year]*n_cells, 
                         [list_metrics]*n_cells, 
                         [list_metrics_kws]*n_cells,
                         [no_data]*n_cells))
    
    # print(poolInput[0])
    
    # Open a pool of processes
    with Pool(num_process) as p:

        print('Starting parallel computation within calculate_classMetrics_withPool.')

        # Sending to compute in parallel
        df_cellMetrics = p.starmap(calculate_classMetrics, poolInput)
        
        # Freeing memory
        del poolInput

        df_cellMetrics = pd.concat(df_cellMetrics, ignore_index=True)

        print('Finished parallel computation within calculate_classMetrics_withPool.')

    return df_cellMetrics


"""
MAIN
"""

def main():
    
    # --------------------------- User inputs --------------------------------
    
    # Asking the user to specifiy the list of years to work on. 
    print('')
    print('------------------------------------------------------------------')
    print('Please specify the range of years to work on \nand the number of processes.')
    print('Input years in the format (separated with space): year1 year2')
    print('')
    server_order = input("Start and end years: \n")
    nprocess = input("Number of Pool processes: \n")
    print('------------------------------------------------------------------')
    print('')
    
    # Casting input as list of years
    list_years = [int(i) for i in server_order.split()]
    list_years = list(range(list_years[0], list_years[1] + 1))
    
    # Casting number of processes as integer
    nprocess = int(nprocess)
    
    # The clipped version of the grid
    grid_fp = '../../data/shapes/cerrado_grid_50km_clipped.shp'
    
    # A test MBlulc file
    # lulc_fp = 'D:/carlota/projects/data/MBlanduse_c06/mapbiomas-brazil-collection-60-cerrado-{}.tif'
    lulc_fp = '../../data/raw/MBlulc/mapbiomas-brazil-collection-60-cerrado-{}.tif'
    
    # Output filepath
    output_fp = '../../data/processed/landscape_metrics/MBlulc_landscapeCellMetrics_grid50km_{}.csv'
    
    # The set of landscape-level metric we want to calculate
    cell_metrics = ['total_area', 'proportion_of_landscape',
                    'number_of_patches', 'patch_density', 'largest_patch_index',
                    'total_edge', 'edge_density', 'landscape_shape_index',
                    'area_mn', 'area_md', 'area_cv']
    cell_metrics_kws = {'patch_density': {'percent': False, 'hectares': True}}
    
    # Threshold to include cells in grid
    area_thresh = 75.0
    
    # The no data value to assign to the raster file
    no_data = 0
    
    
    # ----------------------- Processing data --------------------------------
    
    # Reading and formatting the grid
    grid = read_grid_rasterCRS(grid_fp, lulc_fp, list_years[0], area_thresh)
    
    # grid = grid.head(9)
    
    
    # Working on one cell at a time (to change by parallel)
    for year in list_years:
        
        print('-------------------------------------------------------------')
        print('Working on year {}.'.format(year))
        
        t.tic()
    
        # Calculating the metrics per grid cell
        grid_metrics = calculate_classMetrics_withPool(grid, lulc_fp, year, 
                                                       cell_metrics, cell_metrics_kws, 
                                                       no_data, nprocess)
        
        # Saving dataframe
        print('Saving landscape metrics...\n')
        grid_metrics.to_csv(output_fp.format(year), index = False)
        
        t.toc()
        # break
        
    
    
    
    
if __name__ == '__main__':
    main()