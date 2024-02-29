# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 13:43:50 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculates the area of the polygons of a GeoDataFrame. 

The polygon's area is calculated in a lengthy manner but which assures that the
area will be correct regardles of the shape of the polygons or their validity.
To calculate the area, this script counts the number of MapBiomas Fogo pixels
that fall within each polygon. It transforms number of pixels to units in km2.
The MapBiomas Fogo pixels are the ones used to create the polygons in the first 
instance.

Input files: those found in directories
    ../../data/processed/MBfogo_c10_cerrado_polygons_b150_CM/ (SHP and CSV)
    ../../data/raw/MBfogo_ba/ (TIF)
    
Output files: the CSV files are overwritten
    ../../data/processed/MBfogo_c10_cerrado_polygons_b150_CM/ (CSV)


Functions in this script:
    * calculateArea_km2 - calculates the areas of each polygon in SHP and 
                          writes it as a new column to the corresponding CSV
"""

"""
Imports
"""

import os
import pandas as pd

from pytictoc import TicToc
t = TicToc()

# My modules
import sys
sys.path.append('../utils/')

import fileList
import MB_tools
from assignRasterValue_toPolygon import assignRasterValues_toPolygon
from MBfogo_matchingFileNames import matchingFileNames_year


"""
Function
"""
def calculateArea_km2(gdf_fp, csv_fp, raster_fp):
    """
    Function that calculates the areas of a dataframe of MapBiomas fire polygons by
    counting the number of pixels inside the polygon of value 1 (indicating pixel is burned)
    against 0, no data or not burned

    Parameters
    ----------
    gdf_fp : string
        Filepath to SHP file containing the polygons.
    csv_fp : string
        Filepath to CSV containing each polygon's information.
    raster_fp : string
        Filepath to raw MapBiomas Fogo burned area data.

    Returns
    -------
    None. This function adds the area in km2 as a column ot the CSV file and
    overwrites it.

    """

    
    # 1. Using assignRaster to get the count of 0s and 1s for each polygon.
    # Do this in parallel, as it takes time
    
    print('Counting pixels per polygon')
    t.tic()
    gdf = assignRasterValues_toPolygon(gdf_fp, raster_fp, 0, 12, [0], num_process=35)
    t.toc('Pixels counted in')
    
    
    
    # 2. Formatting informaitn and correcting column names, etc.
    # The function assignRasterValues_toPolygon returns the count of pixels 
    # within polygon for every different possible value of the pixels.
    # In this case, values are 1 to 12, because they indicate the month. Hence,
    # we now have to sum the counts of all the month columns into one,
    # and then erase the extra columns.
    
    # Getting the names of the columns containing the counts per month
    count_columns = list(gdf.columns)[-12:]

    # Calculating total count
    gdf['n_pixels'] = gdf[count_columns].sum(axis=1)
    
    # Erasing extra-columns
    gdf = gdf.drop(columns = count_columns)
    
    
    
    # 3. Converting number of pixels to area (km2)
    # Multiplying pixel counts by pixel area (MB_tools.py)
    # First, select those columns that do not contain the number of pixels - 
    # to pass this information to the function that performs pixel to area transformation
    id_cols = list(gdf.loc[:, gdf.columns != 'n_pixels'].columns)
    
    print('Calculating area.')
    t.tic()
    gdf = MB_tools.pixel_to_area_km2(gdf, id_cols)
    t.toc('Area calculated in')
    
    gdf = gdf.rename(columns={'n_pixels' : 'area'})
    
    

    # 4. Appending this information as anew column to CSV data
    # Convert GeoDataFrame to DataFrame
    gdf = pd.DataFrame(gdf.drop(columns='geometry'))
    
    # Merge this information to the csv table as a column, based on ID
    # First, read CSV file
    gdf_csv = pd.read_csv(csv_fp)
    
    # Merging dataframes
    print('Merging shapefile and csv.')
    gdf_csv = gdf_csv.merge(gdf[['ID', 'area']], on = 'ID')
    
    # Save CSV file overwriting the previous one
    print('Writing CSV file.')
    gdf_csv.to_csv(csv_fp, index=False)
    


"""
MAIN
"""
    
def main():
    
    
    #----------------- User inputs -------------------------------------------
    
    # UPDATE (23/10/2023): MapBiomas collection to work with
    collection = 2
    
    # Filepaths to polygons and their associated data 
    # Directory with input shapefiles and csv, for which we want to calculate 
    # the area
    in_dir = '../../data/processed/MBfogo_c{}0_cerrado_polygons_b150_CM/'.format(collection)
    in_dir = '../../../_emilio/polygon_data_separated/'
    # Shapefiles
    list_shp_fp = fileList.fileList(in_dir, in_extensions='.shp')
    # CSV files
    list_csv_fp = fileList.fileList(in_dir, in_extensions='.csv')


    # Raster files that we use to calculate the area
    raster_dir = '../../data/raw/MBfogo_ba/'
    # raster_dir = 'D:/carlota/projects/data/MBfogo_c20/'
    list_raster_fp = fileList.fileList(raster_dir, in_extensions='.tif')
    # list_raster_fp = [raster_dir + 'MBfogo_c10_cerrado_2000.tif']
    

    
    #-------------- Area calculation -----------------------------------------

    # Preparing the list of filepaths 
    list_fp = matchingFileNames_year([list_shp_fp, list_csv_fp, list_raster_fp])
    
    # Calculating area perpolygon in each set of files
    for shp_fp, csv_fp, raster_fp in list_fp:
        
        t.tic()
        print('Working on file {}'.format(os.path.basename(shp_fp)))
        print(csv_fp)
        print(raster_fp)
        calculateArea_km2(shp_fp, csv_fp, raster_fp)
        t.toc('Area for file {} produced in'.format(os.path.basename(shp_fp)))
        

if __name__ == '__main__':
    main()