# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 10:43:10 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculates a series of landscape-level metrics as a means to 
quantify landscape's spatial heterogeneity. It does so using the package 
"pylandstats", which implements in Python the metrics described in FRAGSTATS.

It calculates the landscpape metrics for each cell and year in the period.


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
def read_grid_rasterCRS(grid_fp, raster_fp, year, thresh = 75.0, include_cells = False, exclude_cells = False):
    
    # Get the raster's CRS
    with rasterio.open(raster_fp.format(year)) as src:
        raster_crs = src.crs
        
    
    # Read the grid
    grid = gpd.read_file(grid_fp)
    
    # Transform to raster's CRS
    grid = grid.to_crs(raster_crs)
    
    if include_cells == False:
        # Filter only cell with x% of the area in Cerrado
        grid = grid[grid['%area_crop'] >= thresh]

    elif include_cells != False:
        # Selecting only those cells whose percentage inside the Cerrado is above a
        # certain threshold, or the cells that are included, and excluding if necessary
        grid = grid[(grid['%area_crop'] >= thresh) | grid['id'].isin(include_cells)].reset_index(drop = True)
        
        if exclude_cells != False:
            grid = grid[~grid['id'].isin(exclude_cells)].reset_index(drop = True)
    
    # # Subsetting relevant column
    # grid = grid[['id']]
    
    return(grid)


def calculate_landscapeMetrics(cell_id, polygon, raster_fp, year, list_metrics, list_metrics_kws, no_data):
    
    # with open("../test_before.txt", "a") as myfile:
    #     myfile.write(str(cell_id)+'\n')
    

    with rasterio.open(raster_fp.format(year)) as raster:
        
        cell_array, cell_transform = mask(dataset = raster, shapes=[polygon], crop=True, nodata = no_data)
        
        # Reducing the number of dimensions (as there is an empty third dimension)
        cell_array = np.squeeze(cell_array)
        

    # Creating the Landscape object
    landsc = pls.Landscape(cell_array, res = (30, 30), nodata = no_data, 
                           transform = cell_transform)

    # Calculate the metrics
    # try:
    metrics_df = landsc.compute_landscape_metrics_df(metrics = list_metrics,
                                                     metrics_kws = list_metrics_kws)
    # except:
    #     with open("../error_in.txt", "a") as myfile:
    #         myfile.write(str(cell_id)+'\n')
    #     raise(Exception('Error with cell {}'.format(cell_id)))
    
    # Including column with year and cell_id
    metrics_df.insert(loc = 0, column = 'year', value = year)
    metrics_df.insert(loc = 0, column = 'id', value = cell_id)
    
    # with open("../test_after.txt", "a") as myfile:
    #     myfile.write(str(cell_id)+'\n')
    
    return metrics_df



def calculate_landscapeMetrics_withPool(grid, raster_fp, year, list_metrics, list_metrics_kws, no_data, num_process = 10):
    
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

        print('Starting parallel computation within calculate_landscapeMetrics_withPool.')

        # Sending to compute in parallel
        df_landscapeMetrics = p.starmap(calculate_landscapeMetrics, poolInput)
        
        # Freeing memory
        del poolInput

        df_landscapeMetrics = pd.concat(df_landscapeMetrics, ignore_index=True)

        print('Finished parallel computation within calculate_landscapeMetrics_withPool.')

    return df_landscapeMetrics


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
    
    grid_id = '30km'
    
    # Casting input as list of years
    list_years = [int(i) for i in server_order.split()]
    list_years = list(range(list_years[0], list_years[1] + 1))
    
    # Casting number of processes as integer
    nprocess = int(nprocess)
    
    # The clipped version of the grid
    # UPDATE (20/10/2023): the cropped version instead! 
    grid_fp = '../../data/shapes/cerrado_grid_{}_cropped.shp'.format(grid_id)
    
    # A test MBlulc file
    # UPDATE (20/10/2023): the MB lulc collection to work with
    lul_collection = 8
    # lulc_fp = 'D:/carlota/projects/data/MBlanduse_c06/mapbiomas-brazil-collection-60-cerrado-{}.tif'
    lulc_fp = '../../data/raw/MBlulc/mapbiomas-brazil-collection-{}0-{}.tif'.format(lul_collection, '{}')
    
    # Output filepath
    output_fp = '../../data/processed/landscape_metrics/MBlulc_landscapeMetrics_grid{}_{}.csv'.format(grid_id, '{}')
    
    # The set of landscape-level metric we want to calculate
    landscape_metrics = ['total_area','number_of_patches', 'patch_density', 'total_edge', 
                          'edge_density', 'landscape_shape_index', 
                          'area_mn', 'area_cv']#, 
                          #conditional_entropy', 'contagion'] # 'entropy', 
    landscape_metrics_kws = {'patch_density': {'percent': False, 'hectares': True}}
    
    # Threshold to include cells in grid
    area_thresh = 0.0 #75.0
    
    # # Cells to include/exclude
    # include_cells = [567, 568]
    # exclude_cells = [1451]
    
    # The no data value to assign to the raster file
    no_data = 0
    
    
    # ----------------------- Processing data --------------------------------
    
    # Reading and formatting the grid
    grid = read_grid_rasterCRS(grid_fp, lulc_fp, list_years[0], area_thresh)#,
                               #include_cells = include_cells, exclude_cells = exclude_cells)
    
    # print(grid)
    
    # Working on one cell at a time (to change by parallel)
    for year in list_years:
        
        print('-------------------------------------------------------------')
        print('Working on year {}.'.format(year))
        
        t.tic()
    
        # Calculating the metrics per grid cell
        grid_metrics = calculate_landscapeMetrics_withPool(grid, lulc_fp, year, 
                                                           landscape_metrics, landscape_metrics_kws, 
                                                           no_data, nprocess)
        
        # Saving dataframe
        print('Saving landscape metrics...\n')
        grid_metrics.to_csv(output_fp.format(year), index = False)
        
        t.toc()
        
    
    
    
    
if __name__ == '__main__':
    main()