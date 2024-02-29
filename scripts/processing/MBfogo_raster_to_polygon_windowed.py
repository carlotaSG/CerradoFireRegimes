# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 10:52:40 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that takes a raster file where pixels values correspond to different 
categories, and clusters adjacent pixels of equal values of the raster into a 
single polygon. Pixels are considered adjacent using a four- or eight-neighbours
criteria, input by the user.

It performs this action over the raster file by windows - small squared polygons -, 
so that it can be parallelised, but also because the function rasterio.features.shapes 
works well in small and simple matrices.

The list of raster files over which to operate is read from a file ./server_orders/MBfogo_raster_to_polygon{}.csv
where {} can be empty or be _node<number>. This allows to divide the work among nodes
and to speify the files over which we wish the script to work on (aka, the years
for which we want to polygonise the MaBiomas [FOGO] data).



"""

"""
Imports
"""
import rasterio
from rasterio.features import shapes
from rasterio.mask import mask

import geopandas as gpd
import pandas as pd

# from shapely import geometry, ops

from multiprocessing import Pool

import numpy as np

from pytictoc import TicToc
t = TicToc()


# Setting directory
import os
# Get full path of current script
abspath = os.path.abspath(__file__)
# Extract directory from the path
dname = os.path.dirname(abspath)
# Change to script directory
os.chdir(dname)


# My modules
import sys
sys.path.append('../utils/')

import fileList
import MBfogo_extract_time_reference as extract_time


# UPDATE 02/10/2023
# I am getting the warning "FutureWarning: pandas.Int64Index is deprecated
# and will be removed from pandas in a future version.Use pandas.Index with 
# the appropriate dtype instead.
# This error seems to be located in /soge-home/users/scat8298/.conda/envs/geo_env/lib/python3.9/site-packages/geopandas/io/file.py:299:
# But I am not keen on updating geopandas. Hence, I am ignoring this specific type
# of warnings.
import warnings

def fxn():
    warnings.warn("deprecated", FutureWarning)



"""
Functions
"""

# TODO: Once I know this works (it runs fine for all years of MBF data), erase this part
# Because polygons are fixed using PyQGIS.
# def resolving_polygon_issues(gdf, cellID):
    
#     # First, fix geometries - we apply .buffer(0) to each of the polygons in gdf
#     # HINT: buffer(0) seems to work more effectively than unary_union(geom)
#     # gdf['geometry'] = gdf['geometry'].apply(lambda geo: ops.unary_union(geo))
#     # # gdf['geometry'] = gdf['geometry'].apply(lambda geo: geo.buffer(0))
    
#     # Then, make sure all geometries are valid
#     invalid = gdf.loc[ ~gdf['geometry'].is_valid, :]
    
#     # If there are any invalid geometries, write to log.
#     if invalid.shape[0] != 0:
#         with open('../outputs/MBfogo_raster_to_polygon_log.txt', 'a') as file:
#             file.write('\n Invalid geometries in cell {} of grid 08 after unary_union.'.format(cellID))
            
    
#     # Finally, checking if any polygons cross one another.
#     # In principle, no polygons should cross one another, not even those of the 
#     # same month.
#     # overlapping = np.sum(gdf['geometry'].apply(lambda x: gdf['geometry'].crosses(x)).values.astype(int))
    
#     # if overlapping > 0 :
#     #     with open('../outputs/MBfogo_raster_to_polygon_log.txt', 'a') as file:
#     #         file.write('\n Overlapping geometries in cell {} of grid 08 after validity check.'.format(cellID))
    
#     # HINT 999: Writing to file in parent function - to be checked by a 
#     # separate script because it takes too long (around 25 minutes per cell, 
#     # out 372 cells or so)
    
#     return gdf
    
    
    
def polygonise(data, month, aff_transform, crs, filterValue = None, connectivity = 8):
    """
    Function that gets a matrix and polygonises it using a certain connectivity.
    Then, it keeps only those polygons with a value in filterValue.

    Parameters
    ----------
    data : np.array
        Matrix of values to be polygonised. It only contains two values, either 0
        or the month that we are working on. Because we are polygonising one month at a time.
    month : integer
        The pixel value that we are going to polygonise.
    aff_transform : rasterio.Affine()
        The affine transformation corresponding to the matrix data, which is given
        to the resulting polygons.
    crs : string
        The Coordinate Reference System of the original matrix data. It is assigned to
        created GeoDataFrame.
    filterValue : list of integers, optional
        The list of pixel values that we wish to keep. If None, all values are 
        kept. The default is None.
    connectivity : integer, optional
        Either 8 or 4, correspond to the neighbouring criteria (number of pixels
        that are considered adjacent). The default is 8.

    Returns
    -------
    polygons : GeoDataFrame
        The list of polygons resulting from polygonising the matrix data. It has 
        the CRS of the original raster and polygon data.contains a column
        called raster_val containing the value of the pixel's clustered together - 
        which is the same as month.

    """
    
    # TODO (2021/12/20): review this. Perhaps, erase. mask is not necessary, and month
    # does not need to be passed. Because the matrix data already only contains the 
    # pixels with value = month and the rest of pixels equal to 0.
    # Create mask (as in https://www.programcreek.com/python/example/118984/rasterio.mask)
    # mask = data == month
    mask = None
    
    # UPDATE 02/10/2023 See update in imports
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fxn()
    
    # Polygonize the matrix
    polygons = shapes(data, mask = mask, transform = aff_transform, connectivity = connectivity)
    
    # Cast it into a GeoJSON-like output
    polygons = list(({'geometry': g, 'properties': {'raster_val': int(v)}} for i, (g, v) in enumerate(polygons)))
    
    # Convert to GeoDataFrame
    polygons = gpd.GeoDataFrame.from_features(polygons)
    
    # Assign raster CRS
    polygons = polygons.set_crs(crs)
    
    # If filterValue is passed as a vector, then we keep only these 'raster_val' values
    if filterValue != None:
        polygons = polygons[polygons['raster_val'].isin(filterValue)]
        
    return polygons


    
def raster_to_polygon(raster_fp, cell_gdf, out_fp, grid_id, year, filterValue = None, connectivity = 8):
    """
    Function that gets a filepath to a raster and a GeoDataFrame with a unique polygon.
    It read the part of the raster inside the polygon and then it polygonises it.
    FInally, it saves it to out_fp.

    Parameters
    ----------
    raster_fp : string
        Filepath to raster file to be parallelised.
    cell_gdf : GeoDataFrame
        Contains a unique polygon used as bounding box to read the raster image.
        cell_gdf CRS must be the same as the raster's CRS.
    out_fp : string
        Filepath to shapefile where to write output polygons.
    grid_id : string
        Indicates the column of grid that identifies each cell uniquely.
    year : integer
        The year of data over which we are working.
    num_process : integer
        The number of workers to use in parallelization.
    filterValue : list of integers, optional
        The list of pixel values that we wish to keep. If None, all values are 
        kept. The default is None.
    connectivity : integer, optional
        Either 8 or 4, correspond to the neighbouring criteria (number of pixels
        that are considered adjacent). The default is 8.

    Returns
    -------
    None.

    """
    
    # Reading raster values found inside cell in grid
    with rasterio.open(raster_fp) as src:
        
        # Read matrix and affine transformation
        raster_data, aff_transform = mask(src, [cell_gdf['geometry']], crop = True)
        
        # Keeping raster CRS in memory to later assign to the created GeoDataFrame
        raster_crs = src.crs

    # If all pixels in cell are 0 (no pixel burned), then write output csv indicating
    # that cell is empty. This helps keeping a log for later, if there are cells missing
    # in shapfile
    if np.sum(raster_data) == 0:
        with open('../logs/MBfogo_raster_to_polygon_log.csv', 'a') as file:
            file.write('\n Empty raster in cell {}.'.format(cell_gdf[grid_id]))
    

    
    
    # We work on the raster only if it is not empty. Otherwise, we jump on to the next one
    if np.sum(raster_data) != 0:
        
        # TODO: (21/09/2021) I am still having the same problem of some pixels
        # not being mapped even using smaller cells (with way smaller cells it does map 
        # more pixels, but some pixels are still left out.)
        # Hence, I am going to try to polygonise by month to see what happens
        
        # List to accumulate the polygons for the different months
        polygons = []
        
        # We loop over the different months:
        for month in range(1, 12 + 1):
            # Subset raster_data to this month
            month_raster_data = np.where(raster_data == month, raster_data, 0)
            
            # Polygonise this subset
            month_polygons = polygonise(month_raster_data, month, aff_transform, raster_crs, [month], connectivity)
            
            # month_polygons is a GeoDataFrame, hence, we accumulate it in a list if it 
            # is not empty, to later on concatenate all of the non-empty months
            if len(month_polygons) > 0:
                polygons.append(month_polygons)
                
        # Converting list of GeoDataFrames to GeoDataFrame
        polygons = gpd.GeoDataFrame( pd.concat(polygons, ignore_index = True))
        
        # Polygonising raster
        #polygons = polygonise(raster_data, aff_transform, raster_crs, filterValue, connectivity)
        
        # Freeing memory
        del raster_data   
        
        # # If the appropriate folder to save this file does not exist, create it
        # folder = os.path.dirname(out_fp).format(year)
        # if not os.path.isdir(folder):
        #     os.makedirs(folder)
        
        # Save polygons to file
        polygons.to_file(folder + '/' + os.path.basename(out_fp).format(grid_id, cell_gdf[grid_id], year))
    
        # Freeing memory
        del polygons

    
    

def raster_to_polygon_windowed(raster_fp, out_fp, grid, grid_id, year, num_process, filterValue = None, connectivity = 8):
    """
    Function that gets a filepath to a raster file to be polygonise and a grid of
    polygons over which to polygonise one geometry at a time (by parts/windows).
    
    It implements this process parallelising over the geometries in grid.

    Parameters
    ----------
    raster_fp : string
        Filepath to raster file to be parallelised.
    out_fp : string
        Filepath to shapefile where to write output polygons.
    grid : GeoDataFrame
        The list of of cells in grid over which to polygonise the raster, one by one.
        Its CRS must already be the same as the raster's CRS.
    grid_id : string
        Indicates the column of grid that identifies each cell uniquely.
    year : integer
        The year of data over which we are working.
    num_process : integer
        The number of workers to use in parallelization.
    filterValue : list of integers, optional
        The list of pixel values that we wish to keep. If None, all values are 
        kept. The default is None.
    connectivity : integer, optional
        Either 8 or 4, correspond to the neighbouring criteria (number of pixels
        that are considered adjacent). The default is 8.

    Returns
    -------
    None.

    """    

    # Prepare input to polygonise in parallel using map
    poolInput = [(raster_fp, grid.iloc[i, :], out_fp, grid_id, year, filterValue, connectivity) for i in range(len(grid.geometry))]

    
    # polygonise raster by parts (cells in grid) and working in parallel
    with Pool(num_process) as p:

        print('Starting parallel computation within raster_to_polygon_windowed.')
        
        # Sending to polygonise in parallel
        # Return list of tuples: GeoDataFrames per cell and cell identifiers
        p.starmap(raster_to_polygon, poolInput)
        
    # Freeing some memory (erasing input to parallelisation)
    del poolInput
    

"""
MAIN
"""
if __name__ == "__main__":
    
    
    # UPDATE 02/10/2023 See update in imports
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fxn()
    
    #---------------------- User inputs -------------------------------
    # Asking the user whether to read server orders for a specific node
    server_order = input("Input '_node<>' to select a specific set of orders.\nFor no-selection, just press enter.\n")
    
    # CSV with the list of raster files to work with in this script
    csv_fp = '../server_orders/MBfogo_raster_to_polygon{}.csv'.format(server_order)   
    
    # # Output directory
    # out_dir = '../data/processed/MBfogo_c10_cerrado_polygons_{}.shp'
    # Generating output filepaths
    out_fp = '../../data/polygons_tocheck/{}/MBfogo_c20_{}{}_{}.shp'
    
    # List of raster values to keep as polygons
    raster_values = list(range(1, 12 + 1))
    
    # Type of connectivity to be used
    connectivity = 8
    
    # Number of processes to be used in parallelization
    num_process = 38
    
    # Read grid file of 0.2 degrees covering the extent of MapBiomas fogo
    # collection 1.0 over the Cerrado. 
    # HINT: The grid's CRS must be the same as the raster's CRS already!!
    # grid_fp = '../../data/shapes/cerrado_grid_02deg_cropped_labelled.shp'
    # grid_id = 'cID02'
    
    # UPDATE 02/10/2023
    # Read grid file of 0.2 degrees covering the extent of a convex hull
    # polygon over the Cerrado
    grid_fp = '../../data/shapes_ch/cerrado_grid_02deg_ch_cropped_labelled.shp'
    grid_id = 'cID02'
    grid = gpd.read_file(grid_fp)
    
    #---------------------- The programme -----------------------------
    
    # Reading list of raster filepaths
    list_raster_fp = fileList.readListFiles_fromCSV(csv_fp)
    
    # Polygonise each raster file in the list
    for in_fp in list_raster_fp:
        print("Working on file {}".format(in_fp))#os.path.basename(in_fp)))
        
        # Extract year of data from filepath
        year = extract_time.extract_year(in_fp)
        
        # If the appropriate folder to save this file does not exist, create it
        folder = os.path.dirname(out_fp).format(year)
        if not os.path.isdir(folder):
            os.makedirs(folder)

        # ---------------- Polygonise -----------------------------------------
        # Polygonise by window: result are files per cell with the polygons checked for validity
        print('Starting polygonising.')
        t.tic()
        raster_to_polygon_windowed(in_fp, out_fp, grid, grid_id, year, num_process, filterValue = raster_values, connectivity = connectivity)
        t.toc('Finished polygonising year {} in'.format(year))


        
    

