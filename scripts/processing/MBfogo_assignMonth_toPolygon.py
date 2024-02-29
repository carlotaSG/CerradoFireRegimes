# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 11:20:56 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk


Script that gets a GeoDataFrame of polygons and assigns a month label to each
polygon. Returns a GeoDataFrame with the format
    year | month

The month assigned is the minimum raster value of the original MapBiomas Fogo
data that is found within the polygon. This task is completed by the script
../utils/minimumRasterValue_toPolygon.py

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

import fileList
from minimumRasterValue_toPolygon import minimumRasterValues_toPolygon


"""
Functions
"""

def get_fp(list_years, shp_dir, rast_dir):
    """
    Function that return a list of tuples containing the corresponding pairs of
    (shapefile, rasterfile) for each year in list_years. It does so by first reading
    all files in each directory - either shp or tif.

    Parameters
    ----------
    list_years : list of integers
        List of years for which we want to retrieve the filepaths.
    shp_dir : string
        Path to directory containing the shapefiles of interest.
    rast_dir : string
        Path to directory containing the rasterfiles of interest.

    Returns
    -------
    shape_raster_fp : list of tuples of three positions: year, shapefilepath, raster filepath
        Contains the pairs of filepaths to the shapefile and rasterfile corresponding
        to each year in the list, along with the respective year.

    """
    
    # Reading list of all shapefiles in directory
    shapefiles_fp = fileList.fileList(shp_dir, in_extensions='.shp')

    
    # Reading list of rasterfiles
    rasterfiles_fp = fileList.fileList(rast_dir, in_extensions='.tif')
    
    
    # Selecting only those shapefiles/rasterfiles corresponding to the list of years
    shapefiles_fp = fileList.selecting_fp_withCondition(shapefiles_fp, list_years)  
    rasterfiles_fp = fileList.selecting_fp_withCondition(rasterfiles_fp, list_years)
    
    # Creating tuples of the pairs
    shape_raster_fp = list(zip(list_years, shapefiles_fp, rasterfiles_fp))
    
    return shape_raster_fp
    


"""
MAIN
"""

def main():
    
    
    # Printing WARNING message that this script overwrites original files
    print('-------------------------------------------------')
    print('WARNING: This script overwrites the input shapefiles in:')
    print('../../data/processed/MBfogo_c10_cerrado_polygons_b150_CM/MBfogo_c10_cerrado_polygons_{}.{}')
    print('')
    user_decision = input('Are you OK with this? (y/n)\n')
    print('-------------------------------------------------')
    
    if user_decision == 'n':
        sys.exit(1)
    
    # ------------------------- User inputs -----------------------------------
    
    # Asking the user whether to read server orders for a specific node
    server_order = input("Input '_node<>' to select a specific set of orders.\nFor no-selection, just press enter.\n")
    
    # CSV with the list of raster files to work with in this script
    csv_fp = '../server_orders/MBfogo_assignMonth_toPolygon{}.csv'.format(server_order) 
    
    
    # Reading list of files to work on (one after the other)
    list_years = fileList.readListFiles_fromCSV(csv_fp)
    
    # MapBiomas collection to work with
    collection = 2
    
    
    # Declaring directories
    # Containing the shapefiles
    shapefiles_dir = '../../data/processed/MBfogo_c{}0_cerrado_polygons_b150_CM/'.format(collection)
    
    # Containing the raster files
    rasterfiles_dir = '../../data/raw/MBfogo_ba/'
    
    
    # ---------------------- Getting filepaths --------------------------------
    
    # Getting the pairs of filepaths needed to complete the task. 
    # That is, the shapefiles containing the polygons and the rasterfiles from 
    # which to get the month label
    # shape_raster_fp = get_fp(list_years, shapefiles_dir, rasterfiles_dir)
    
    
    # print(shape_raster_fp)
    #------------------- The core part -----------------------------------
    
    # Now, we loop over each shapefile-rasterfile pair
    for year in list_years:
        
        print('----------------------------------------------------------------')
        print('Working on year {}'.format(year))
        
        # Shapefile
        shapefile_fp = shapefiles_dir + 'MBfogo_c{}0_cerrado_polygons_{}.shp'.format(collection, year)
        
        # Raster file
        rasterfile_fp = rasterfiles_dir + 'MBfogo_c{}0_cerrado_{}.tif'.format(collection, year)
        
        print(shapefile_fp)
        print(rasterfile_fp)
        
        # Get the minimum raster value within each polygon
        t.tic()
        gdf_min = minimumRasterValues_toPolygon(shapefile_fp, rasterfile_fp, nodata = 0, num_process=29)
        t.toc('Month label was assigned in ')
        
        # Format the GeoDataFrame
        gdf_min = gdf_min[['year', 'min_val', 'geometry']]
        
        # Rename column containing the minimum month to month
        gdf_min = gdf_min.rename(columns = {'min_val': 'month'})
        
        # Order by month
        gdf_min = gdf_min.sort_values(by = ['year','month'])
        

        # Saving GeoDataFrame to Shapefile
        t.tic()
        gdf_min.to_file(shapefile_fp)
        t.toc('SHP file was saved in ')
        
        print('')
        
        
    
if __name__ == '__main__':
    main()
