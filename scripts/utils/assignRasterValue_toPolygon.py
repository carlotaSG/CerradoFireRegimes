# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 15:13:40 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that, given a GeoDataFrame of polygons and a raster file, counts the number 
of raster pixels within each polygon and pixel category. Return a GeoDataFrame
with the initial information and the number of pixels per pixel category.

At the moment, the pixels categories should be integers.

This file can also be imported as a module and contains the following functions:
    * countRasterValues_inPolygon - counts number of pixels per pixel category
                                    within a polygon
    * countRasterValues_inListPolygons - counts number of pixels per pixel category
                                         given a list of polygons
    * countRasterValues_inListPolygons_withPool - prepares the input to run the
                                    function countRasterValues_inListPolygons
                                    in parallel. Parallelising over the list of 
                                    polygons
    * assignRasterValues_toPolygon - Function that calls the above functions
                                    for a certain GeoDataFrame and reaster file,
                                    returning the GeoDataFrame with the counts of 
                                    pixels per pixel category.
    * createGDF - NOT USED TODO: move to another script?
    
"""

"""
Imports
"""
# Setting directory
import os
import sys
# # Get full path of current script
# abspath = os.path.abspath(__file__)
# # Extract directory from the path
# dname = os.path.dirname(abspath)
# # Change to script directory
# os.chdir(dname)

#import rasterstats as rstats
import numpy as np
import rasterio
from rasterio.mask import mask
import geopandas as gpd
import pandas as pd

from pytictoc import TicToc
t = TicToc()

import fileList
from multiprocessing import Pool


"""
Functions
"""
def countRasterValues_inPolygon(polygon, raster_fp, min_val, max_val, nodata):
    """
    Function that receives a raster filepath and a polygon and returns the number
    of pixels within the polygon by pixel value. Raster values must be discrete.
    Returns the counts per value in range [min_val, max_val].
    
    Rater values must hence be numerical.
    
    It will actually return the counts for all values in the range [min_val, max_val],
    even if some values in range do not appear in raster.
    
    IMPORTANT: the polygon's CRS must already be the same as the raster's CRS.

    Parameters
    ----------
    polygon : Shapely geometry
        Polygon geometry for which the raster pixels are going to be counted.
    raster_fp : String
        Filepath to raster file from which the raster values counts are to be obtained.
    min_val : integer/float
        Minimum possible value of raster pixels.
    max_val : integer/float
        Maximum possible value of raster pixels.

    Returns
    -------
    counts : TYPE
        DESCRIPTION.

    """

    # Open raster file in 'r' mode
    with rasterio.open(raster_fp) as raster:
        
        # If nodata is None, then the raster's nodata is used
        if nodata == None:
            # Read only the raster information within the polygon
            out_img, out_transform = mask(dataset = raster, shapes=[polygon], crop=True)
        
        # Otherwise, use the nodata provided by user
        elif nodata != None:
            out_img, out_transform = mask(dataset = raster, shapes = [polygon], crop = True, nodata = nodata)
        
    # Creating bins for each one of the classes (numbers) between min_val and max_val
    # Each class/number is the in the middle of the bin, hence adding/subtracting 0.5
    binZ = np.arange(min_val - 0.5, max_val + 1 + 0.5, 1).tolist()
    
    # Counting how many pixels are in each bin
    counts, bins = np.histogram(out_img, bins=binZ)
    
    return counts


def createGDF(polygon_fp, min_category, max_category):
    """
    Creates the GeoDataFrame where the raster value counts per polygon are going
    to be stored. Returned dataframe contains all of the information of the input
    GeoDataFrame (polygon_fp), and empty columns - one for each value in range
    [min_category, max_category].
    
    Note that it will create a column for each value in this range, regardless of
    whether that value exists in raster.

    Parameters
    ----------
    polygon_fp : string
        Filepath to the GeoDataFrame that contains all polygons for which to count
        raster values of each pixel class, and extra column with additional information.
    min_category : integer/float
        Prefereably integer, minimum value the raster pixels can take. Used to 
        create the sequence of columns of raster values.
    max_category : integer/float
        Prefereably integer, maximum value the raster pixels can take. Used to 
        create the sequence of columns of raster values.

    Returns
    -------
    gdf : GeoDataFrame
        GeoDataFrame contianing all information in input GeoDataFrame (polygon_fp),
        as well as empty columns for each value in the range [min_val, max_val].
    category_columns : List of strings
        List of the column names corresponding to each value in the range [min_val, max_val].

    """

    # Read the GeoDataFrame
    polygon_gdf = gpd.read_file(polygon_fp)
    
    # Create list of column names
    polygon_columns = polygon_gdf.columns.tolist()
    category_columns = [str(cat) for cat in list(range(min_category, max_category + 1))]
    columns = polygon_columns + category_columns
    
    # Create empty GeoDataFrame with these columns
    gdf = gpd.GeoDataFrame(columns= columns)
    
    # Fill the polygon_columns with the data from polygon_gdf
    gdf[polygon_columns] = polygon_gdf[polygon_columns].copy()
    
    return gdf, category_columns
    

def countRasterValues_inListPolygons(input_vars):
    """
    Creates a GeoDataFrame with the counts of the number of pixels per pixel 
    value for a list of polygons (or a unique polygon). Raster values must be discrete.
    Returns the counts per value in range [min_val, max_val].
    
    It will actually return the counts for all values in the range [min_val, max_val],
    even if some values in range do not appear in raster.
    
    IMPORTANT: the polygon's CRS must already be the same as the raster's CRS.

    Parameters
    ----------
    gdf_geometries : List of geometries or GeoPandas GeoSeries.
        List of polygon geometries for which the raster values are going to be counted.
    raster_fp : String
        Filepath to raster file from which the raster values counts are to be obtained.
    min_val : integer/float
        Minimum possible value of raster pixels.
    max_val : integer/float
        Maximum possible value of raster pixels.

    Returns
    -------
    df_categories : Pandas DataFrame
        DataFrame with the number of pixels within polygon per pixel value category.

    """
    # Unpacking input
    gdf_geometries, raster_fp, min_val, max_val, nodata = input_vars
    
    # If gdf_geometries is not a list of geometries, cast into one
    if isinstance(gdf_geometries, gpd.geoseries.GeoSeries) == False:
        gdf_geometries = [gdf_geometries] 
    
    
    # Create empty dataframe one column per value in the range [min_category, max_category]
    #    List of the range of categories as column names
    category_columns = [str(cat) for cat in list(range(min_val, max_val + 1))]
    
    #    Create empty GeoDataFrame with these columns
    # df_categories = gpd.GeoDataFrame(columns= category_columns)
    
    # Create list to accumulate dataframe
    df_categories = []

    # For each polygon in list gdf_geometries
    for polygon in gdf_geometries:
        
        # Count number of raster values within polygon
        counts = countRasterValues_inPolygon(polygon, raster_fp, min_val, max_val, nodata)

        # Cast into dataframe
        polygon_df_categories = pd.DataFrame([counts], columns=category_columns)
        
        # Append it to the cumulative list of DataFrames
        df_categories.append(polygon_df_categories)
        
        # Free memory
        del polygon_df_categories
    
    # Cast list of dataframes to dataframe
    df_categories = pd.concat(df_categories, ignore_index = True)

    
    return df_categories

def countRasterValues_inListPolygons_withPool(gdf_geometries, raster_fp, min_val, max_val, nodata, num_process):
    """
    Prepares the input to run the function countRasterValues_inListPolygons in parallel
    with a number of processes num_process. It parallelises over the list of
    polygon geometries.

    Parameters
    ----------
    gdf_geometries : list oh Shapely geometries
        List of geometries for which to count the number of raster values per 
        category.
    raster_fp : String
        Filepath to raster file.
    min_val : integer/float
        Minimum possible value of raster pixels.
    max_val : integer/float
        Maximum possible value of raster pixels.
    nodata : integer/float
        Value to be used as nodata assignment when reading raster. If None, then
        the raster's nodata value will be used instead.
    num_process : integer
        Number of processes to create use in the parallelisation.

    Returns
    -------
    df_rasterValues : Pandas Dataframe
        DataFrame with the number of pixels within polygon per pixel value category
        for every geometry in gdf_geometries.

    """

    # Preparing input for Pool
    n_polygons = len(gdf_geometries)
    poolInput = list(zip(gdf_geometries, 
                         [raster_fp]*n_polygons, 
                         [min_val]*n_polygons, 
                         [max_val]*n_polygons,
                         [nodata]*n_polygons))
    
    # Open a pool of processes
    with Pool(num_process) as p:

        print('Starting parallel computation within countRasterValues_inListPolygon_withPool.')

        # Sending to compute in parallel
        df_rasterValues = p.map(countRasterValues_inListPolygons, poolInput)
        
        # Helping to free memory (erasing input to parallelisation)
        del poolInput

        df_rasterValues = pd.concat(df_rasterValues, ignore_index=True)

        print('Finished parallel computation within countRasterValues_inListPolygon_withPool.')

    return df_rasterValues

def assignRasterValues_toPolygon(polygon_fp, raster_fp, min_val, max_val, nodata = None, non_categories = None, num_process=None):
    """
    Function that goes through every step from receiving a filepath for a 
    GeoDataFrame of geometries and for a raster file, to extracint the pixel 
    counts per pixel value category for every geometry (polygon) in the polygon_fp
    GeoDataFrame.

    Parameters
    ----------
    polygon_fp : string
        Fielpath to the GeoDataFrame containing the list of polygons for which
        to retrieve the pixel counts.
    raster_fp : string
        Filepath to the raster file from which to count the number of pixels.
    min_val : integer/float
        Minimum possible value of raster pixels.
    max_val : integer/float
        Maximum possible value of raster pixels.
    nodata : integer/float, optional
        Value to be assigned as 'nodata' when reading and masking the raster. If not 
        provided, script will use the metadata information. If this information is
        None too, script will inform the user on command line and will terminate the
        programme.
    non_categories : list of integers, optional
        List of those values within the range [min_val, max_val] that do not
        exist in the raster. The default is None.
    num_process : integer, optional
        Number of processes with which to parallelise. I None, then the polygons
        are processed recurrently. The default is None.

    Returns
    -------
    gdf : GeoDataFrame
        GeoDataFrame with the same information in polygon_fp plus the count of
        the number of pixels per pixel value category for each value in range
        [min_val, max_val], except those values passed in non_categories.

    """
    
    # 1. Read the polygon GeoDataFrame
    gdf = gpd.read_file(polygon_fp)
    
    # Its CRS
    gdf_initial_crs = gdf.crs
    
    
    with rasterio.open(raster_fp) as raster:
        # If user has not provided a nodata value, and the raster's metadata is 
        # None too, programme informs user and terminates
        if nodata == None:
            if raster.nodata == None:
                print("WARNING:")
                print("--------------------------------------------------------------")
                print("nodata value of raster file: \n {}".format(raster_fp))
                print("is None and user has not provided an alternative nodata value.")
                print("--------------------------------------------------------------")
                print("Terminating programme. Nodata value is required")
                sys.exit(1)
        
        # Adjusting the CRS of the GeoDataFrame to that of the raster file
        raster_crs = raster.crs
        gdf = gdf.to_crs(raster_crs)
    

    # 2. Now we pass the geometry column to a function that will retrieve the counts per pixel
    # category for each polygon (each row in dataframe). It will return a dataframe 
    # with the pixel categories as columns.
    
    # Consecutively
    if num_process == None:
        df_rasterValues = countRasterValues_inListPolygons((gdf['geometry'], raster_fp, min_val, max_val, nodata))

    # In parallel
    elif num_process != None:
        df_rasterValues = countRasterValues_inListPolygons_withPool(gdf['geometry'], raster_fp, min_val, max_val, nodata, num_process)  
    
    
    # 3. Then, we append horizontally gdf_rasterValues to gdf with the CRS of the raster
    gdf = gpd.GeoDataFrame(pd.concat([gdf, df_rasterValues], axis=1), crs = raster_crs)
    
    # 4. Convert back to the original crs
    gdf = gdf.to_crs(gdf_initial_crs)
    

    # 5. Now we have the full GeoDataFrame with the counts of pixels per polygon for this shapefile.

    # Erasing those columns corresponding to values that are not categories of 
    # pixel values if a list of values is passed
    if non_categories != None:
        if isinstance(non_categories[0], str) == False:
            non_categories = [str(element) for element in non_categories]
        # Dropping corresponding categories
        gdf = gdf.drop(columns=non_categories)
    
    
    # Fill NaNs values with 0s because Fiona does not deal with NaNs (I think)
    gdf = gdf.fillna(0)
    
    
    return gdf


"""
MAIN
"""
def main():
    
    #------------------------ User input --------------------------------
        
    # Directory with shapefiles
    shapefiles_dir = '../outputs/MODIS_BA_polygons/annual_pixels_countInstances/'
       
    # Directory with raster files
    rasterfiles_dir = '../data/mapbiomas_LU/cerrado/'
    
    # List of years that we are interested in (to filter out the files we are not interested in)
    years = list(range(2001,2019 + 1))       # TODO: modify this line to send to different nodes
    years = [str(year) for year in years]
    
    # Range of cateogries in the rasterfiles (including non-categories if they fall in range)
    min_cat = 0 # Because this is the category that will be assigned to pixels with no data
    max_cat = 41 
    
    #range_cat = list(range(min_cat, max_cat + 1))
    noncat = [0, 6,7,8,16,17,28,34,35,37,38,40]
    
    
    
    #-------------------- Reading lists of filepaths --------------------
    
    # Reading list of shapefiles
    shapefiles_fp = fileList.fileList(shapefiles_dir, in_extensions='.shp')

    
    # Reading list of rasterfiles
    rasterfiles_fp = fileList.fileList(rasterfiles_dir, in_extensions='.tif')
    
    
    # Selecting only those shapefiles/rasterfiles corresponding to the list of years
    shapefiles_fp = fileList.selecting_fp_withCondition(shapefiles_fp, years)  
    rasterfiles_fp = fileList.selecting_fp_withCondition(rasterfiles_fp, years)
    
    # Creating tuples of the pairs
    shape_raster_fp = list(zip(shapefiles_fp, rasterfiles_fp))
   # print(shape_raster_fp)
    
    
    #------------------- The core part -----------------------------------
    
    # Now, we loop over each shapefile-rasterfile pair
    # Initialising a year iterator
    year_iterator = 0
    for shapefile_fp, rasterfile_fp in shape_raster_fp:
        
        # Working on year
        year = years[year_iterator]
        
        print('Working on \n {}'.format(shapefile_fp))
        
        # We create the corresponding GeoDataFrame with the shapefile_fp information and the 
        # count of pixels per polygon and land-use category
        t.tic()
        gdf_LU = assignRasterValues_toPolygon(shapefile_fp, rasterfile_fp, min_cat, max_cat, noncat, num_process=20)
        t.toc('LU was assigned in ')
        

        # Saving GeoDataFrame to Shapefile
        t.tic()
        gdf_LU.to_file('../outputs/MODIS_BA_polygons/annual_pixels_countInstances/with_LU/MODIS_BA_{}.shp'.format(year))
        t.toc('SHP file was saved in ')
        
 
        # Increasing year_iterator
        year_iterator += 1
        

if __name__ == '__main__':
    main()
