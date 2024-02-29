# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 15:13:40 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that, given a GeoDataFrame of polygons and a raster file, returns the
certain summary raster value within each polygon. Returns a GeoDataFrame with the initial 
information and the summary pixel value (mean, min, max).


This file can also be imported as a module and contains the following functions:
    * quantityRasterValues_inPolygon - counts number of pixels per pixel category
                                    within a polygon
    * quantityRasterValues_inListPolygons - counts number of pixels per pixel category
                                         given a list of polygons
    * quantityRasterValues_inListPolygons_withPool - prepares the input to run the
                                    function countRasterValues_inListPolygons
                                    in parallel. Parallelising over the list of 
                                    polygons
    * quantityRasterValues_toPolygon - Function that calls the above functions
                                    for a certain GeoDataFrame and reaster file,
                                    returning the GeoDataFrame with the counts of 
                                    pixels per pixel category.
    
"""

"""
Imports
"""
# Setting directory
import sys

import numpy as np
import rasterio
from rasterio.mask import mask
import geopandas as gpd
import pandas as pd
from multiprocessing import Pool

from pytictoc import TicToc
t = TicToc()




"""
Functions
"""
def quantityRasterValues_inPolygon(polygon, raster_fp, nodata = None, quantity = 'mean'):
    """
    Function that receives a raster filepath and a polygon and returns the 
    quantity pixel value within the polygon. By default, quantity is the average value.
    
    Raster values must hence be numerical.
    
    IMPORTANT: the polygon's CRS must already be the same as the raster's CRS.

    Parameters
    ----------
    polygon : Shapely geometry
        Polygon geometry for which the minimum raster values are to be obtained.
    raster_fp : String
        Filepath to raster file from which the minimum raster values are to be obtained.
    nodata: numeric, default = None
        The value to be used as no data when reading the raster. If not indicated,
        then the raster's value is used instead.
    quantity: String
        The summary amount to calculate for the polygon.

    Returns
    -------
    val : Numeric
        The quantity raster value inside the polygon.

    """

    # Open raster file in 'r' mode
    with rasterio.open(raster_fp) as raster:
        
        # If nodata is None, then the raster's nodata is used
        if nodata == None:
            # Read only the raster information within the polygon
            out_img, out_transform = mask(dataset = raster, shapes=[polygon], crop=True, filled = True)
            
            # Retrieve nodata value from raster
            nodata = raster.nodata
        
        # Otherwise, use the nodata provided by user
        elif nodata != None:
            out_img, out_transform = mask(dataset = raster, shapes = [polygon], crop = True, nodata = nodata, filled = True)
    
    
    # Flatten the array
    out_img = out_img.flatten()
    
    # Remove nodata values
    out_img = out_img[out_img != nodata]
    
    # Calculating the value to return
    if quantity == 'mean':
        # Calculate mean
        val = np.mean(out_img)
        
    elif quantity == 'min':
        val = np.min(out_img)
        
    elif quantity == 'max':
        val = np.max(out_img)
    
    
    return val
    

def quantityRasterValues_inListPolygons(input_vars):
    """
    Creates a GeoDataFrame with the quantity pixel value for a list of polygons 
    (or a unique polygon).
    
    IMPORTANT: the polygon's CRS must already be the same as the raster's CRS.

    Parameters
    ----------
    gdf_geometries : List of geometries or GeoPandas GeoSeries.
        List of polygon geometries for which the minimum raster values are to be obtained.
    raster_fp : String
        Filepath to raster file from which the minimum raster values are to be obtained.
    nodata: numeric, default = None
        The value to be used as no data when reading the raster. If not indicated,
        nodata should be None and then the raster's value is used instead.
    quantity : String
        Keyword indicating the quantity to be calculated (mean, min, max)


    Returns
    -------
    df_val : Pandas DataFrame
        DataFrame with the quantity pixel value for each polygon.

    """
    # Unpacking input
    gdf_geometries, raster_fp, nodata, quantity = input_vars
    
    # If gdf_geometries is not a list of geometries, cast into one
    if isinstance(gdf_geometries, gpd.geoseries.GeoSeries) == False:
        gdf_geometries = [gdf_geometries] 
    
    
    # Create empty list to accumulate DataFrames
    df_val = []
    

    # For each polygon in list gdf_geometries
    for polygon in gdf_geometries:
        
        # Minimum raster value within polygon
        val = quantityRasterValues_inPolygon(polygon, raster_fp, nodata, quantity)

        # Cast into dataframe
        polygon_val = pd.DataFrame([val], columns=['val'])
        
        # Append it to the cumulative GeoDataFrame
        df_val.append(polygon_val)
        
        # Free memory
        del polygon_val
    
    # Convert list of DataFrames to unique DataFrame
    df_val = pd.concat(df_val, ignore_index = True)

    
    return df_val

def quantityRasterValues_inListPolygons_withPool(gdf_geometries, raster_fp, nodata, quantity, num_process):
    """
    Prepares the input to run the function minimumRasterValues_inListPolygons in parallel
    with a number of processes num_process. It parallelises over the list of
    polygon geometries.

    Parameters
    ----------
    gdf_geometries : list oh Shapely geometries
        List of geometries for which to count the number of raster values per 
        category.
    raster_fp : String
        Filepath to raster file.
    nodata : integer/float
        Value to be used as nodata assignment when reading raster. If None, then
        the raster's nodata value will be used instead.
    quantity : String
        Keyword indicating the quantity to be calculated (mean, min, max)
    num_process : integer
        Number of processes to create use in the parallelisation.

    Returns
    -------
    df_rasterValues : Pandas Dataframe
        DataFrame with the quantity pixel value within polygon for every geometry 
        in gdf_geometries.

    """

    # Preparing input for Pool
    n_polygons = len(gdf_geometries)
    poolInput = list(zip(gdf_geometries, 
                         [raster_fp]*n_polygons, 
                         [nodata]*n_polygons,
                         [quantity]*n_polygons))
    
    # Open a pool of processes
    with Pool(num_process) as p:

        print('Starting parallel computation within quantityRasterValues_inListPolygons_withPool.')

        # Sending to compute in parallel
        df_rasterValues = p.map(quantityRasterValues_inListPolygons, poolInput)
        
        # Helping to free memory (erasing input to parallelisation)
        del poolInput

        df_rasterValues = pd.concat(df_rasterValues, ignore_index=True)

        print('Finished parallel computation within quantityRasterValues_inListPolygons_withPool.')

    return df_rasterValues



def quantityRasterValues_toPolygon(polygon_fp, raster_fp, nodata = None, quantity = None, num_process = None):
    """
    Function that goes through every step from receiving a filepath for a 
    GeoDataFrame of geometries and for a raster file, to extracing the quantity 
    pixel value for every geometry (polygon) in the polygon_fp GeoDataFrame.

    Parameters
    ----------
    polygon_fp : string
        Fielpath to the GeoDataFrame containing the list of polygons for which
        to retrieve the pixel counts.
    raster_fp : string
        Filepath to the raster file from which to count the number of pixels.
    nodata : integer/float, optional
        Value to be assigned as 'nodata' when reading and masking the raster. If not 
        provided, script will use the metadata information. If this information is
        None too, script will inform the user on command line and will terminate the
        programme.
    quantity : String
        Keyword indicating the quantity to be calculated (mean, min, max)
    num_process : integer, optional
        Number of processes with which to parallelise. I None, then the polygons
        are processed recurrently. The default is None.

    Returns
    -------
    gdf : GeoDataFrame
        GeoDataFrame with the same information in polygon_fp plus the quantity
        pixel value.

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
        df_rasterValues = quantityRasterValues_inListPolygons((gdf['geometry'], raster_fp, nodata, quantity))

    # In parallel
    elif num_process != None:
        df_rasterValues = quantityRasterValues_inListPolygons_withPool(gdf['geometry'], raster_fp, nodata, quantity, num_process)  
    
    
    # 3. Then, we append horizontally gdf_rasterValues to gdf with the CRS of the raster
    gdf = gpd.GeoDataFrame(pd.concat([gdf, df_rasterValues], axis=1), crs = raster_crs)
    
    # 4. Convert back to the original crs
    gdf = gdf.to_crs(gdf_initial_crs)
    

    # 5. Now we have the full GeoDataFrame with the minimum pixel values per 
    # polygon for this shapefile.
    
    
    return gdf


"""
MAIN
"""
def main():
    
    # NOTE: Code commented out below comes from a much older version of this script,
    # probably from assignRasterValue_toPolygon. Left for legacy if script needs
    # to be as main.
    print('Nothing to show...')
    
    # #------------------------ User input --------------------------------
        
    # # Directory with shapefiles
    # shapefiles_dir = '../../data/temp/fixed/2000_CM/'
       
    # # Directory with raster files
    # rasterfiles_dir = '../../data/raw/MBfogo_ba/'
    
    # # List of years that we are interested in (to filter out the files we are not interested in)
    # years = [2000] 
    # years = [str(year) for year in years]
    
    
    
    # #-------------------- Reading lists of filepaths --------------------
    
    # # Reading list of shapefiles
    # shapefiles_fp = fileList.fileList(shapefiles_dir, in_extensions='.shp')

    
    # # Reading list of rasterfiles
    # rasterfiles_fp = fileList.fileList(rasterfiles_dir, in_extensions='.tif')
    
    
    # # Selecting only those shapefiles/rasterfiles corresponding to the list of years
    # shapefiles_fp = fileList.selecting_fp_withCondition(shapefiles_fp, years)  
    # rasterfiles_fp = fileList.selecting_fp_withCondition(rasterfiles_fp, years)
    
    # # Creating tuples of the pairs
    # shape_raster_fp = list(zip(shapefiles_fp, rasterfiles_fp))
    # print(shape_raster_fp)
    
    
    # #------------------- The core part -----------------------------------
    
    # # Now, we loop over each shapefile-rasterfile pair
    # # Initialising a year iterator
    # year_iterator = 0
    # for shapefile_fp, rasterfile_fp in shape_raster_fp:
        
    #     # Working on year
    #     year = years[year_iterator]
        
    #     print('Working on \n {}'.format(shapefile_fp))
        
    #     # We create the corresponding GeoDataFrame with the shapefile_fp information and the 
    #     # minimum pixel polygon and land-use category
    #     t.tic()
    #     gdf_min = minimumRasterValues_toPolygon(shapefile_fp, rasterfile_fp, nodata = 0, num_process=10)
    #     t.toc('LU was assigned in ')
        
    #     print(gdf_min)
    #     # Saving GeoDataFrame to Shapefile
    #     t.tic()
    #     gdf_min.to_file('../../data/temp/with_month/{}.shp'.format(year))
    #     t.toc('SHP file was saved in ')
        
 
    #     # Increasing year_iterator
    #     year_iterator += 1
        

if __name__ == '__main__':
    main()
