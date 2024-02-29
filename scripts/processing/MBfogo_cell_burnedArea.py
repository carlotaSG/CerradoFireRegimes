# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 13:43:50 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculates the burned area within a grid's cell. 

The burned area is calculated in a lengthy manner but which assures that the
area will be correct regardless of the shape of the polygons or their validity.
To calculate the area, this script counts the number of MapBiomas Fogo pixels
that fall within each polygon cell. It transforms number of pixels to units in km2.
The MapBiomas Fogo pixels are the ones used to create the polygons in the first 
instance.

Input files: those found in directories
    ../../data/shapes/cerrado_grid_<grid>_cropped.shp (SHP and CSV)
    ../../data/raw/MBfogo_ba/MBfogo_c<collection>0_cerrado_<year>.tif (TIF)
    
Output files: the CSV files are overwritten
    ../../data/processed/summary_tables/<nyears>_year_periods/burnedArea_raster_grid<grid_id>_period<p>.csv (CSV)


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

# import fileList
import MB_tools
from assignRasterValue_toPolygon import assignRasterValues_toPolygon
# from MBfogo_matchingFileNames import matchingFileNames_year


"""
Function
"""
def calculateArea_km2(gdf_fp, raster_fp, nprocess):
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
    gdf = assignRasterValues_toPolygon(gdf_fp, raster_fp, 0, 12, [0], num_process = nprocess)
    t.toc('Pixels counted in')
    
    
    
    # 2. Formatting information and correcting column names, etc.
    # The function assignRasterValues_toPolygon returns the count of pixels 
    # within polygon for every different possible value of the pixels.
    # In this case, values are 1 to 12, because they indicate the month. Hence,
    # we now have to sum the counts of all the month columns into one to have total 
    # area, as well as per month
    
    # Getting the names of the columns containing the counts per month
    count_columns = list(gdf.columns)[-12:]

    # Calculating total count
    gdf['n_pixels'] = gdf[count_columns].sum(axis=1)
    

    
    # 3. Converting number of pixels to area (km2)
    # Multiplying pixel counts by pixel area (MB_tools.py)
    # First, select those columns that do not contain the number of pixels - 
    # to pass this information to the function that performs pixel to area transformation
    id_cols = list(gdf.loc[:, ~gdf.columns.isin(['n_pixels'] + count_columns)].columns)
    # id_cols = list(gdf.loc[:, gdf.columns != 'n_pixels'].columns)
    
    print('Calculating area.')
    t.tic()
    gdf = MB_tools.pixel_to_area_km2(gdf, id_cols)
    t.toc('Area calculated in')
    
    gdf = gdf.rename(columns={'n_pixels' : 'area_T'})
    
    

    # 4.Converting this information to a dataframe
    # Convert GeoDataFrame to DataFrame
    gdf = pd.DataFrame(gdf.drop(columns='geometry'))
    
    # 5. Return the information
    return(gdf)    


"""
MAIN
"""
    
def main():
    
    
    #----------------- User inputs -------------------------------------------
    
    # Filepaths to polygons and their associated data 
    # Directory with input shapefiles and csv, for which we want to calculate 
    # the area
    in_dir = '../../data/shapes/'
    # Shapefiles
    list_grids = ['30km']

    # UPDATE (24/10/2023): inlcude MapBiomas fogo collection option
    collection = 2

    # Raster files that we use to calculate the area
    raster_fp = '../../data/raw/MBfogo_ba/MBfogo_c{}0_cerrado_{}.tif'.format(collection, '{}')
    # list_raster_fp = fileList.fileList(raster_dir, in_extensions='.tif')
    # list_raster_fp = [raster_dir + 'MBfogo_c10_cerrado_2000.tif']
    
    # Groups of years - script will generate a summary file per group
    group_years = [(1985, 1993), (1994, 2002), (2003, 2011), (2012, 2020)]
    
    # Number of years in period
    nyears = 9
    
    output_fp = '../../data/processed/summary_tables/{}_year_periods/burnedArea_raster_grid{}_period{}.csv'
    
    annual_output_fp = '../../data/processed/summary_tables/annual/burnedArea_raster_grid{}_year{}.csv'
    
    # Number of processes to create
    n_process = 35
    
    #-------------- Area calculation -----------------------------------------
    
    # First, looping over each grid we want to calculate the area for
    for grid_id in list_grids:
        
        grid_fp = in_dir + 'cerrado_grid_{}_cropped.shp'.format(grid_id)
        
        # Counter to loop over the periods
        counter = 1
        
        # Looping over groups of years
        # We work on a group of years at a time
        for (start_year, end_year) in group_years:
            
            period = '({}, {})'.format(start_year, end_year)
            
            print('Working on group ' + period)
            
            # Creating empty list to accumulate the burned area for each cell in each year and month
            df_burnedArea_period = []
            
        
            # For each year in group, calculate the burned area inside each cell
            for year in range(start_year, end_year + 1):
                print('')
                print('Working on year {}'.format(year))
    
                # The raster file for this year
                raster_fp_year = raster_fp.format(year)
                
                t.tic()
                df_burnedArea_year = calculateArea_km2(grid_fp, raster_fp_year, n_process)
                t.toc('Area for file {} produced in'.format(os.path.basename(raster_fp_year)))
                
                # -----------------------------------------------------------------
                # UPDATE 17/11/2023: I need the annual information per cell!!
                
                # Create column with year id
                df_burnedArea_year['year'] = year
                # Casting columns in order
                df_burnedArea_year = df_burnedArea_year[['id', 'year','area_T'] + [str(m) for m in range(1, 12+1)]]
                # Saving to file
                df_burnedArea_year.to_csv(annual_output_fp.format(grid_id, year), index = False)
                # Dropping year column
                df_burnedArea_year = df_burnedArea_year.drop(columns = 'year')
                # -----------------------------------------------------------------
                
                # Accumulating this information for the period
                df_burnedArea_period.append(df_burnedArea_year)
                
                # Freeing memory
                del df_burnedArea_year
                
            # Converting list of dataframes to dataframes
            df_burnedArea_period = pd.concat(df_burnedArea_period, ignore_index = True)
            

            # Finally This dataframe contains duplicated cells (one per year)
            # Summarise
            df_burnedArea_period = df_burnedArea_period.groupby(by = 'id').sum().reset_index()

            # Now we add a column to each indicating the period
            df_burnedArea_period['period'] = period
            
            # Adding column to indicate the grid 
            df_burnedArea_period['grid_id'] = grid_id

            
            # Reshuffling columns
            df_burnedArea_period = df_burnedArea_period[['id','period','area_T'] + [str(m) for m in range(1, 12+1)]]
            
            print('Saving data for all years in period...\n')
            
            # Saving dataframes
            df_burnedArea_period.to_csv(output_fp.format(nyears, grid_id, counter), index = False)
            
            counter += 1
            
            
        

if __name__ == '__main__':
    main()