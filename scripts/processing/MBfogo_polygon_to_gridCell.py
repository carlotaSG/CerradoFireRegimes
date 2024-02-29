# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 12:17:37 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that assigns each fire polygon in a shapefile to its corresponding cell
in a grid. The output is a dataframe of lists of polygons per cell. The criteria 
used to assigna polygon to a specific cell is: a polygon is assigned to the cell
where most of the polygon's area falls.

Input files:
    '../../data/shapes/cerrado_grid_<grid_size_units>_cropped.shp'
    '../../data/processed/MBfogo_c10_cerrado_polygons_b150_CM/MBfogo_c10_cerrado_polygons_<year>.<>'
    '../../data/raw/MBfogo_ba/MBfogo_c10_cerrado_<year>.tif'
    

Output files
    ../../data/processed/MBfogo_c10_cerrado_polygons_b150_CM/MBfogo_c10_cerrado_polygons_<year>.csv
    pol_id | year | month | cell_assign
    
where pol_id is the unique identifier of the fire polygon, cell_assign
is the cell to which it has been assigned, and year and month correspond to the 
year and month in which the fire was mapped.

Instead of using geopandas.geometry.area, this script import assignRasterValue_toPolygon,
which counts the number of pixels per pixel category inside each polygon TODO

*** The original script was archived on 17/03/2023. The only changes between that 
version and the current one have to do with tyding up: cleaning the code, commenting
and documenting.


This script has functions:
    * assign_polygon_to_cell - assigns each polygon to the grid cell where most
                               of the polygon's area falls in

"""

"""
Imports
"""
import geopandas as gpd
import pandas as pd
import numpy as np
from multiprocessing import Pool
from pytictoc import TicToc
t = TicToc()

import sys
sys.path.append('../utils')

import fileList
import assignRasterValue_toPolygon
import parallel_prep as pprep


"""
Functions
"""           
        
def assign_polygon_to_cell(slices, polygon_fp, raster_fp, grid, grid_size_units):
    """
    Function that assigns each polygon in polygon_fp to a cell of the grid.
    It does so only for the polygons in polygon_fp whose index is within slices.
    
    It considers that a polygon belongs to the cell where most of its area falls.
    For each polygon, it gets the cells intersecting the polygon. If more than
    one cell is intersecting, it caluclates the intersection area between the 
    polygon and each cell (using assignRasterValues_toPolygon) and selects the 
    cell with the largest intersection. If there is only one cell intersecting, 
    that is the cell returned.
    
    
    Parameters
    ----------
    slices : tuple
        First element indicates the first row in polygons_fp to be read, second element
        indicates the last row to be read.
    polygon_fp : string
        Filepath to the shapefile containing the list of polygons to be assigned
        to a cell.
    raster_fp : string
        Filepath to the raster file contaiing the MBfogo data, to count how many 
        polygon's pixels fall in each cell and calculate the intersection area.
    grid : GeoDataFrame
        GeoPandas dataframe with the list of cells in grid.
    grid_size_units : string
        String that specifies the grid size and units used as: <grid_size><units>
        E.g. 50km
        

    Returns
    -------
    list_df : DataFrame
        DataFrame with the list of polygons together with the cell they have 
        been assigned to. Dataframe returned has columns
        'ID' | cellID_{grid_size_units} 

    """
    
    # Read chunk of polygons in polygon_fp file
    polygons = gpd.read_file(polygon_fp, rows = slice(slices[0], slices[1]))
    
    # Convert grid CRS to that of polygon
    grid = grid.to_crs(polygons.crs)
    
    # Create empty list of dataframes to accumulate each polygon's identification
    list_df = []
    
    # Loop over polygons to assign them to cells individually
    for i in range(len(polygons)): 
        
        # Selecting only the polygon geometry and id
        pol_geom = polygons.loc[i, 'geometry']
        pol_id = polygons.loc[i, 'ID']
        
        
        # First, get the list of cells intersecting the polygon
        cell_inters = grid.loc[grid['geometry'].intersects(pol_geom), 'id'].apply(int).reset_index(drop = True)

        
        
        # If there is just one cell intersecting the polygon
        if len(cell_inters) == 1:
            # Then, this cell and polygon are the ones that go to the list of
            # 'cell_assign'
            cell_assign = cell_inters.item()
            
        # If the polygon intersect more than one cell    
        elif len(cell_inters) > 1:
            
            # In this case, we calculate the geometries of each cell intersection
            cell_assign_geom = grid.loc[grid['id'].isin(cell_inters), 'geometry'].intersection(pol_geom).reset_index(drop = True)
            
            # Then, we calculate the area of each intersection
            # For this we use the following function, which counts the number of
            # pixels inside each geometry per pixel category
            pixels_per_geom = assignRasterValue_toPolygon.countRasterValues_inListPolygons(
                (cell_assign_geom, raster_fp, 1, 12, 0))
            
            # This returns a dataframe with the results for each geometry in each row
            # Then, I just sum across all columns to have a list of pixels per geometry
            pixel_per_geom = pixels_per_geom.sum(axis = 1)
            
            # Then, get the index of cell that has the largest area inside 
            # polygon (largest number of pixels)
            pixel_per_geom = pixel_per_geom.idxmax()
            
            # Use it to subset the cell_ID
            cell_assign = cell_inters[pixel_per_geom].item()
            
        # UPADTE (20/10/2023): For MapBiomas collection 2.0 we have generated the 
        # polygons for a big area including the Cerrado and a conisderable buffer. 
        # Hence, there will be polygons that may not be assigned to any cell.
        # In that case, we assign the cell id -99
        elif len(cell_inters) == 0:
            cell_assign = -99

            
        # Create dataframe with this information for Criteria 2
        polygon_cell = pd.DataFrame({'ID' : pol_id,
                                               'cellID_{}'.format(grid_size_units) : [cell_assign]})
        # Append to list of dataframes obtained using Criteria 2
        list_df.append(polygon_cell)
        
    # Freeing memory (2022/01/20)
    del polygons, polygon_cell
        
    # Cast list of dataframes as a unique dataframe
    list_df = pd.concat(list_df, ignore_index = True)

    
    return list_df


"""
MAIN
"""

def main():
    
    # ----------------------- User inputs -------------------------------------
    
    # Asking the user to specifiy the server order file from which to read the 
    # grid size to use and the list of years
    print('')
    print('------------------------------------------------------------------')
    print('Please select the server_order file that you want to work on.')
    print('The server order file MUST contain <grid_size><units> in the first row,')
    print('followed by the list of years we want the script to work on.')
    print('')
    server_order = input("Input '_node<>' to select a specific set of orders.\nFor no-selection, just press enter.\n")
    print('------------------------------------------------------------------')
    print('')
    
    # MapBiomas collection we are working with
    collection = 2
    
    # CSV with the list of raster files to work with in this script
    csv_fp = '../server_orders/MBfogo_polygon_to_gridCell{}.csv'.format(server_order) 
    
    # Reading list of files to work on (one after the other)
    list_years = fileList.readListFiles_fromCSV(csv_fp)

    # Grid size
    grid_size_units = list_years[0]
    
    # Converting years from string to year
    list_years = [int(year) for year in list_years[1:]]
    
    # Number of processes to parallelise over
    num_process = 8
    
    # ------ Classifying polygon by cell for each year in period -------------- 
    
    
    # Grid filepath
    grid_fp = '../../data/shapes/cerrado_grid_{}_cropped.shp'.format(grid_size_units)
    # Reading grid
    grid = gpd.read_file(grid_fp)
    
    
    # For each year, classify polygons
    for year in list_years:
        print('')
        print('Working on year {}'.format(year))
    
        # ---------------------------------------------------------------------
        # Getting paths to input files
    
        # Filepath to files
        files = '../../data/processed/MBfogo_c{}0_cerrado_polygons_b150_CM_fixed/MBfogo_c{}0_cerrado_polygons_{}.{}'.format(collection, collection, '{}', '{}')
        # Filepath to fire polygon shapefile
        polygon_fp =  files.format(year, 'shp')
        
        # Filepath to fire polygon data
        polygon_data_fp = files.format(year, 'csv')
        
        # Filepath to MBfogo raster
        raster_fp = '../../data/raw/MBfogo_ba/MBfogo_c{}0_cerrado_{}.tif'.format(collection, year)
        
        # # The number of processes to parallelise over
        # num_process = 10
    
        # ---------------------------------------------------------------------
        # Assigning each polygon to a cell
    
        # We assign the polygons to a certain cell by parallelising over groups 
        # of polygons. To achieve this, we first calculate the groups of polygons 
        # as slices of indices (row names). Each process in Pool will be in 
        # charge of a certain group of rows of the polygons file.
        
        # Reading the polygon's csv file to get the total number of rows.
        polygon_data = pd.read_csv(polygon_data_fp)
        # Only the total number of rows is needed to be able to produce partitions
        npolygons = int(len(polygon_data.index))
        del polygon_data
        
        # Generate the slices
        slices = pprep.list_slices(npolygons, num_process*4)
        
        
        # Generating input to Pool using the slices
        pool_input = list(zip(slices, [polygon_fp]*len(slices), [raster_fp]*len(slices), 
                              [grid]*len(slices), [grid_size_units]*len(slices)))

    

        t.tic()
        with Pool(num_process) as p: 
            
            print("Starting parallel computation in assign_polygon_to_cell...")
            
            # Assign polygon to a certain cell, return a dataframe ? TODO
            df = p.starmap(assign_polygon_to_cell, pool_input)
            
            print('Finished parallel computation.')
            
            
        # Converting list of dataframes to a single dataframe
        df = pd.concat(df, ignore_index = True)
        
        
        # Then, we read again the csv datafile and add cellID as a new column (2022/01/20)
        # Read file
        polygon_data = pd.read_csv(polygon_data_fp)
        # Merge new column based on 'ID' (the polygons' ID)
        polygon_data = pd.merge(polygon_data, df, on = 'ID', how = 'outer')
        
        # Freeing memory
        del df
        
        print('Writing data to output files.')
        polygon_data.to_csv(
            polygon_data_fp,
            index = False)
    
    t.toc('Polygon-cell matching done in')
        
        
if __name__ == '__main__':
    main()

