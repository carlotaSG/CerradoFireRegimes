# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 12:17:37 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

***Update 2022/01/10: This script only generates data for criteria 2. Code lines
    criteria 1 have been commented out.

***Update 2022/01/20: This script has been modified so that the cell classification
    using criteria 2 is added as a column to the dataframe in ../data/processed/
    MBfogo_c10_cerrado_polygons/MBfogo_c10_cerrado_polygons_{year}.csv, which is
    overwritten.
    
***Update 2022/01/20: Moved function list_slices() to independent script 
   ../../scripts/parallel_prep.py so that it can be used by other scripts

Script that assigns each fire polygon in a shapefile to its corresponding cell
in a grid. The output is a dataframe of lists of polygons per cell. There are two
criteria to assign a polygon to a grid cell: (1) if the polygon simply intersects
- considering a buffer to account for float point precision - the cell polygon; 
(2) if the majority of the polygon's area falls inside the cell. The second 
is more restrictive.

There are two output files, one for each criteria. Each one contains the relation
between each polygon and the grid cell(s) to which it belongs depending on the criteria.

Criteria 1:
    data/processed/grid_{}km/MBfogo_polygon_grid{size}km_cerrado_{year}_criteria01.csv
    pol_id | month | cell_inters
    
Criteria 2:
    data/processed/grid_{}km/MBfogo_polygon_grid{size}km_cerrado_{year}_criteria02.csv
    pol_id | month | cell_assign

where pol_id is the unique identifier of the fire polygon, cell_inters/assign
is the cell to which it has been assigned, and month is the month in which the fire was mapped.

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
import assignRasterValue_toPolygon
import parallel_prep as pprep


"""
Functions
"""
# def list_slices(ntotal, partitions):
#     """
#     Function that, given the total number of rows in file and a number of partitions,
#     returns the slices of the range (0, ntotal), as many slices as the number of partitions.
    
#     All slices have the same length except the last one, which includes the
#     mode (the residual) in its length.   

#     Parameters
#     ----------
#     ntotal : integer
#         The numer of rows in the data (the range of values).
#     partitions : integer
#         The number of partitions in which to divide the rows (range of values).

#     Returns
#     -------
#     slices : list of tuples
#         List of slices of indices (row names) to be used when partitioning the 
#         data for parallel computation.

#     """
    
#     # Number of elements per slice rounded down
#     nelements = int(ntotal/partitions)
    
#     # The mode or residual
#     the_mod = ntotal % partitions
    
#     # The right hand side of the slice
#     right = np.arange(0, nelements * partitions, nelements)
    
#     # The left hand side is the same as the right but starting from the second element
#     left = right[1:len(right)]
#     # And adding the very last element
#     left = np.append(left, left[len(left)-1] + nelements + the_mod)
    
#     # Finally, the slices as a list of tuples
#     slices = list(zip(right, left))
    
#     return slices
            
        
def assign_polygon_to_cell(slices, polygon_fp, raster_fp, grid, grid_size, grid_units):
    """
    Function that assigns each polygon in polygon_fp to a cell of the grid.
    It does so only for the polygons in polygon_fp whose index is within slices.
    
    It does so using two criteria: (1) a polygon belongs to a cell if it touches it;
    (2) a polygon belongs to the cell where most of its area falls.
    
    ***Update 2022/01/10 only criteria 2 used, criteria 1 commented out.
    ***Update 2022/01/20 Changed name of output dataframe columns, which implied
                         adding two new arguments to function.

    Parameters
    ----------
    slices : tuple
        First element indicates the first row in polygons_fp to be read, second element
        indicates the last row to be read.
    polygon_fp : string
        Filepath to the shapefile containing the list of polygons to be assigned
        to a cell.
    raster_fp : string
        Filepath to the raster file contaiing the MBfogo data, to be used by 
        criteria 2 to count how many polygon's pixels fall in each cell.
    grid : GeoDataFrame
        GeoPandas dataframe with the list of cells in grid.
    grid_size : int (2022/01/20)
        The length of the grid cells' side
    grid_units : str (2022/01/20)
        The units of grid_size

    Returns
    -------
    list_df_criteria2 : DataFrame
        DataFrame with the list of polygons together with the cell they have 
        been assigned to according to criteria 2. Dataframe returned has columns
        'ID' | cellID_{grid_size}{grid_units} (Updated 2022/01/20)

    """
    
    # Read chunk of polygons in polygon_fp file
    polygons = gpd.read_file(polygon_fp, rows = slice(slices[0], slices[1]))
    
    # Convert grid CRS to that of polygon
    grid = grid.to_crs(polygons.crs)
    
    # Create empty list of dataframes
    # list_df_criteria1 = []
    list_df_criteria2 = []
    
    # Loop over polygons to assign them to cells individually
    for i in range(len(polygons)): 
        
        # Polygon geometry, id and month
        pol_geom = polygons.loc[i, 'geometry']
        pol_id = polygons.loc[i, 'ID']
        #pol_month = polygons.loc[i, 'month']   2022/01/20 commented out
        
        
        # First, list of cells intersecting the polygon
        cell_inters = grid.loc[grid['geometry'].intersects(pol_geom), 'id'].apply(int).reset_index(drop = True)

        # # 1. Criteria 1:
        # # Then, create dataframe with this information
        # polygon_cell_criteria1 = pd.DataFrame({'pol_id' : pol_id,
        #                                 'month' : pol_month,
        #                                 'cell_inters': cell_inters})
        
        
        # # Append to list of dataframes obtained using Criteria 1
        # list_df_criteria1.append(polygon_cell_criteria1)
        
        
        # 2. Criteria 2:
        
        # If there is just one cell intersecting the polygon
        if len(cell_inters) == 1:
            # Then, this cell and polygon are the ones that go to the list of
            # 'cell_assign'
            cell_assign = cell_inters.item()
        elif len(cell_inters) > 1:
            
            # In this case, we calculate the geometries of each cell intersection
            cell_assign_geom = grid.loc[grid['id'].isin(cell_inters), 'geometry'].intersection(pol_geom).reset_index(drop = True)
            
            # Then, we calculate the area of each intersection
            pixels_per_geom = assignRasterValue_toPolygon.countRasterValues_inListPolygons(
                # gdf_geometries, raster_fp, min_val, max_val, nodata
                (cell_assign_geom, raster_fp, 1, 12, 0))
            
            # This returns a dataframe with the results for each geometry in each row
            # Then, I just sum across all columns to have a list of pixels per geometry
            pixel_per_geom = pixels_per_geom.sum(axis = 1)
            # Index of cell that has has the largest area inside polygon (largest
            # number of pixels)
            pixel_per_geom = pixel_per_geom.idxmax()
            
            # Use it to subset the cell_ID
            cell_assign = cell_inters[pixel_per_geom].item()

            
        # Create dataframe with this information for Criteria 2
        polygon_cell_criteria2 = pd.DataFrame({'ID' : pol_id,
                                               #'month' : pol_month,   2022/01/20 Commented out
                                               'cellID_{}{}'.format(grid_size, grid_units) : [cell_assign]})   # 2022/01/20 cell_assign -> 'cellID_{}{}'.format(grid_size, grid_units)
        # Append to list of dataframes obtained using Criteria 2
        list_df_criteria2.append(polygon_cell_criteria2)
        
    # Freeing memory (2022/01/20)
    del polygons
        
    # Cast list of dataframes as a unique dataframe
    # list_df_criteria1 = pd.concat(list_df_criteria1, ignore_index = True)
    list_df_criteria2 = pd.concat(list_df_criteria2, ignore_index = True)

    
    # return (list_df_criteria1, list_df_criteria2)
    return list_df_criteria2


"""
MAIN
"""

def main():
    
    # ----------------------- User inputs -------------------------------------
    # Grid size
    grid_size = 50
    # Units of grid_size
    grid_units = 'km'
    
    # Period over which to work
    # Start year
    start_year = 2000
    # Final year
    end_year = 2000
    
    # ------ Classifying polygon by cell for each year in period -------------- 
    # Filepath to grid
    grid_fp = '../../data/shapes/cerrado_grid_{}{}_cropped.shp'.format(grid_size, grid_units)
    # Reading grid
    grid = gpd.read_file(grid_fp)
    
    # For each year, classify polygons
    for year in range(start_year, end_year + 1):
        print('')
        print('Working on year {}'.format(year))
    
        # Filepath to files
        files = '../../data/processed/MBfogo_c10_cerrado_polygons_b150_CM/MBfogo_c10_cerrado_polygons_{}.{}'
        # Filepath to fire polygon shapefile
        polygon_fp =  files.format(year, 'shp')
        # polygon_fp = '../../data/tests/fire_polygons/nobuff_consecutiveMonths/MBfogo_c10_cerrado_polygons_subset_{}_nobuff.shp'.format(year)
        
        # Filepath to fire polygon data
        polygon_data_fp = files.format(year, 'csv')
        # polygon_data_fp = '../../data/tests/fire_polygons/nobuff_consecutiveMonths/MBfogo_c10_cerrado_polygons_subset_{}_nobuff.csv'.format(year)
        
        # Filepath to MBfogo raster
        raster_fp = '../../data/raw/MBfogo_ba/MBfogo_c10_cerrado_{}.tif'.format(year)
        
        # The number of processes to parallelise over
        num_process = 10
    
    # ------------------------------------------------------------------------
    # Notes (30/11/2021)
    
    # Now, I could use Dask Geopandas to run operations in parallel using Dask's
    # functionalities. I could also use MultiProcessing with Pool. Not sure which one
    # is faster in reading the shapefile.
    # Apparently, Dask does not read the shapefile. What you do is read the shapefile 
    # in Geopandas and then import it to Dask, which partitions it by rows and 
    # just makes parallelising over it easier.
    
    
    # As well, I wonder if it might just be easier to read over small chunks 
    # of the shapefile, instead of all at once. Not sure what would be faster.
    # That is, I would get the number of rows from reading the length of the 
    # csv  and then I would partition over slices of the shapefile. Reading
    # each bit at a time and doing the operations over this small chunk.
    
    # I think that reading small chunks of a dataframe at a time might be a 
    # better option, because it takes shorter and because I could be reading 
    # different chunks of the shapefile at the same time, hence reducing the
    # time of overall reading.
    
    # ------------------------------------------------------------------------
    
    
        # We assign the polygons to a certain cell by parallelising over groups 
        # of polygons. To achieve this, we first calculate the groups of polygons 
        # as slices of indices (row names). Each process in Pool will be in 
        # charge of a certain group of rows of the polygons file.
        
        # Reading the polygon's csv file to get the total number of rows.
        polygon_data = pd.read_csv(polygon_data_fp)
        npolygons = int(len(polygon_data.index))
        del polygon_data
        
        # We just want the total number of rows to be able to produce partitions of
        # the file
        
        # Generate the slices
        slices = pprep.list_slices(npolygons, num_process*4)
        
        # Generating input to Pool using the slices
        pool_input = list(zip(slices, [polygon_fp]*len(slices), [raster_fp]*len(slices), 
                              [grid]*len(slices), [grid_size]*len(slices), [grid_units]*len(slices)))

    

        t.tic()
        with Pool(num_process) as p: 
            
            print("Starting parallel computation in assign_polygon_to_cell...")
            
            # output = p.starmap(assign_polygon_to_cell, pool_input)
            df_criteria2 = p.starmap(assign_polygon_to_cell, pool_input)
            
            print('Finished parallel computation.')
            
            # It returns a list of tuples with (df_criteria1, df_criteria2)
            # df_criteria1 = [element[0] for element in output]
            # df_criteria2 = [element[1] for element in output]
            
            
            # Converting list of dataframes to a single dataframe
            # df_criteria1 = pd.concat(df_criteria1, ignore_index=True)
            df_criteria2 = pd.concat(df_criteria2, ignore_index = True)
            
            # print(len(df_criteria1))
            # print(len(df_criteria2))
            
            # Then, we read again the csv datafile and add cellID as a new column (2022/01/20)
            # Read file
            polygon_data = pd.read_csv(polygon_data_fp)
            # Merge new column based on 'ID' (the polygons' ID)
            polygon_data = pd.merge(polygon_data, df_criteria2, on = 'ID', how = 'outer')
            
            print('Writing data to output files.')
            # df_criteria1.to_csv('../data/processed/MBfogo_polygon_grid{}km_cerrado_{}_criteria01.csv'.format(grid_size, year), index = False)
            polygon_data.to_csv(
                polygon_data_fp.format(year),
                index = False)
            # 2022/01/20 Changed ../data/processed/grid_{}km/MBfogo_polygon_grid{}km_cerrado_{}_criteria02.csv
            #            by ../data/processed/MBfogo_c10_cerrado_polygons/MBfogo_c10_cerrado_polygons_{}.csv
            
            
            
            # print(pd.concat(g for _, g in df_criteria1.groupby("pol_id") if len(g) > 1))
            # print(pd.concat(g for _, g in df_criteria2.groupby("pol_id") if len(g) > 1))
            
            # Code seems to be working fine and giving reasonable reults!
        
        
        t.toc('Grid-cell matching done in ')
if __name__ == '__main__':
    main()

