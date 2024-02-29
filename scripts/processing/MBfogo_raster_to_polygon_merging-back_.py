# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 10:52:40 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that takes polygons given in different shapefiles each belonging to a cell
in a grid and merges those polygons touching one another - within a certain 
buffer - and with the same polygon value. It does so recurrently using the same
process with various grids of increasing cell size until recovering the whole 
of the Cerrado.

It does so following the grid order indicated in ../aid_inputs/raster_to_polygon_merging_info.csv
which contains the information of the grid on which to work, it's intersections (crosses),
the input folder from which to read the files, and the output folder to where write
the resulting the files. This script writes to a temporarily folder ../temp/
to diminish RAM usage.



"""

"""
Imports
"""
# import rasterio
# from rasterio.features import shapes
# from rasterio.mask import mask

# import sys

import geopandas as gpd
import pandas as pd

# from shapely import geometry, ops

from multiprocessing import Pool

# import numpy as np
from scipy.sparse.csgraph import connected_components

# import matplotlib.pyplot as plt

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
# import MB_tools
import MBfogo_extract_time_reference as extract_time

import warnings
warnings.filterwarnings('ignore')

"""
Functions
""" 
def dissolve_polygons(polygons, cell_id):
    """
    Function that dissolves those polygons in dataframe that are touching one 
    another (merges them).
    Code from https://gis.stackexchange.com/questions/271733/geopandas-dissolve-overlapping-polygons
    
    Merging is done using a buffer of 0.0001

    Parameters
    ----------
    polygons : GeoDataFrame
        List of polygons touching one another to be merged.
    cell_id : string
        Identifier of the cell that we are merging to.

    Returns
    -------
    polygons : GeoDataFrame
        List of polygons resulting from dissolving/merging those that are touching
        one another.

    """

    # Matrix that indicates which polygons area touching one another
    try:
        touching_matrix = polygons['geometry'].apply(lambda x: polygons['geometry'].overlaps(x.buffer(0.0001))).values.astype(int)
        # (24/09/2021) Changing buffer from 0.00028 to 0.0001
    except:
        with open('../logs/dissolve_polygons_log.csv', 'a') as file:
                    file.write('\n index error in cell  {}.'.format(cell_id))
    # TODO: solving problem: Reindexing only valid with uniquely valued Index objects
    # It seems to be related to the function apply over a dataframe that has a repeated index...
    # It can be solved by using .reset_index(drop=True), but first I want to know where this happens
    
    #print(touching_matrix)
    # Reformatting information to get a series of ids
    n, ids = connected_components(touching_matrix)
    
    # Add this information as a column back to the polyogns GeoDataFrame
    polygons['group'] = ids
    #print(ids)
    
    
    
    
    # Then, dissolve by group
    polygons = polygons.dissolve(by = 'group').reset_index(drop=True)
    
    # TODO: trying to solve the lines that are left over
    # First, enlarging geometries to make sure inner lines will be erased.
    # Using join_style = 2 (mitre, pointy ends)
    polygons['geometry'] = polygons['geometry'].buffer(0.0001, join_style = 2)
    # Finally, shrink back geometries again
    polygons['geometry'] = polygons['geometry'].buffer(-0.0001, join_style = 2)
    
    return polygons
        

def merge_polygons(arg_input):
    """
    

    Parameters
    ----------
    arg_input : iterable of inputs
        Intersections to be used to select cnadidate polygons to be merged,
        Filepaths to shapefiles that are to be merged together (they all fall in
        next level grid cell),
        Next grid ID to be used as file code to identify the erged polygons,
        Year in which we are working,
        Folder path to write output polygons to.

    Returns
    -------
    None.

    """
    
    # GeoDataFrame of lines/crosses in cell
    cross = arg_input[0].reset_index(drop=True)
    cross = cross.loc[0, 'geometry']
    #print(cross)
    
    # List of polygon files in cell
    polygons_fp = arg_input[1]
    
    # # ID of cell in grid of smaller size
    # small_grid_ID = arg_input[2]
    #ID of cell in grid of larger size
    large_grid_ID = arg_input[2]
    
    #folder = arg_input[4]
    
    # Year 
    year = arg_input[3]
    
    # Output folder
    output_folder = arg_input[4]
    
    # If there are no polygon files produced for this large_grid_ID cell, we skip this cell
    # Otherwise, we merge the polygons inside
    if len(polygons_fp) == 0:
        with open('../logs/MBfogo_raster_to_polygon_log.csv', 'a') as file:
                file.write('\n polygons_fp empty for cell {}'.format(large_grid_ID))
    if len(polygons_fp) > 0:
        
        # Read files with polygons and concatenate to unique GeoDataFrame
        polygons = []
        for fp in polygons_fp:
            #print(fp)
            polygons.append(gpd.read_file(fp))    
            
        # Inroducing ignore_index = True to avoid duplicate indices
        polygons = gpd.GeoDataFrame( pd.concat(polygons, ignore_index=True) )

        #print(polygons.head())
        
        # TESTING: If there are no intersections, len(polygons_fp) should be one. In that case,
        # there are no intersections and so there are no polygons to dissolve (merge). The only
        # required action is to rewrite the file with the new file naem (corresponding to the next grid)
        # and save to corresponding folder.
        if len(polygons_fp) == 1:
            all_polygons = polygons.copy()
            
            # Writing log message
            with open('../logs/MBfogo_raster_to_polygon_log.csv','a') as file:
                file.write('\n Only one smaller cell for a large cell in cell {}'.format(large_grid_ID))
            
        elif len(polygons_fp) > 1:
            
        
            # List to accumulate the resulting gdfs of polygons per month
            dissolved_polygons = []
            
            # For every month, we subset the corresponding polygons, and then those that
            #are touching the cross
            for month in range(1, 12 + 1):
                polygons_month = polygons[polygons['raster_val'] == month]
                # print(len(polygons_month))
                
                try:
                    # Subsetting polygons touching cross lines
                    polygons_cross_month = polygons_month[polygons_month['geometry'].intersects(cross.buffer(0.00028))]
                    # print(len(polygons_cross_month))
                except:
                    with open('../logs/MBfogo_raster_to_polygon_log.csv', 'a') as file:
                        file.write('\n polygons_cross_month error  {}'.format(large_grid_ID))
                        for i in range(len(polygons_fp)):
                            file.write('\n ' + polygons_fp[i])
                    sys.exit(1)
                # Subsetting polygons not touching cross lines
                polygons_notcross_month = polygons_month[~polygons_month['geometry'].intersects(cross.buffer(0.00028))]
                # print(len(polygons_notcross_month))
                
                # Freeing memory
                del polygons_month
                
                # Merging into one polygon those polygons touching the cross that are 
                # touching one another.
                if len(polygons_cross_month) != 0:
                    polygons_cross_month = dissolve_polygons(polygons_cross_month, large_grid_ID)
                    
                polygons_month = gpd.GeoDataFrame( pd.concat( [polygons_notcross_month, polygons_cross_month], ignore_index = True ))
                
                # print('\n Checking polygons_month {}'.format(month))
                # print(polygons_month.head())
                
                dissolved_polygons.append(polygons_month)
                
                # print(len(dissolved_polygons))
                
                #Freeing memory
                del polygons_month
                
                
                # break
            
            # Concatenating again the gdfs for each month
            all_polygons = gpd.GeoDataFrame( pd.concat(dissolved_polygons, ignore_index = True) )
        
        #print(polygons.head())
        
        # Save to file
        # If the appropriate folder to save this file does not exist, create it
        #folder = '../data/temp/{}/'.format(year)
        folder = output_folder.format(year)
        if not os.path.isdir(folder):
            os.makedirs(folder)
        all_polygons.to_file(folder + 'MBfogo_c20_{}_{}.shp'.format(large_grid_ID, year))
        
        
        #return polygons
    
        # print(type(grid_lines.loc[0, 'geometry'].intersection(grid_lines.loc[3, 'geometry'])))
        # inter = gpd.GeoSeries(grid_lines.loc[0, 'geometry'].intersection(grid_lines.loc[2, 'geometry']))
        # inter.plot()
        # plt.show()
        
        # TODO: erase files in polygon_fp


def merging_polygons_in_grid(initial_grid, input_polygon_folder, output_polygon_folder, 
                             grid_intersections_fp, grid_id, next_grid_id, year, num_processes):
    """
    Function that merges all polygons in a certain grid that belong to the same
    month and different cells that are touching one another. We do this in parallel
    operating over mutliple cells of a larger polygon at the same time. That is,
    operating over groups of cells in current cells, groups are defined by cells 
    in the next level grid.

    Parameters
    ----------
    initial_grid : GeoDataFrame
        The current grid in which the ucrrent polygons are distributed in.
    input_polygon_folder : string
        Filepath to folder containing the polygons in current grid, to be merged.
    output_polygon_folder : string
        Filepath to folder where to write the polygons merged and distributed
        according to next grid.
    grid_intersections_fp : string
        Filepath to GeoDataFrame containing the crosses or cell intersections of
        current grid. These are used to select those polygons that are candidates
        to be merged.
    grid_id : string
        Identifier of current grid.
    next_grid_id : string
        Identifier of next grid, by which to merge the polygons in current grid.
    year : int
        Year of data we are working on.
    num_process : int
        Number of processes to use in parallelization.

    Returns
    -------
    None.

    """
    
    print('Working on merging polygons')
    
    # Getting list of unique cell IDs in grid
    cell_id = initial_grid[next_grid_id].unique()
    # HINT: since we are getting the list of larger-cell grid ids from the initial
    # grid, if there are larger cells that contain no smaller cells, these will not
    # be listed to be chosen in the pool.
    
    # Get list of files with polygons in grid 
    folder = input_polygon_folder.format(year)
    files = fileList.fileList(folder, in_extensions = 'shp')
    
    # Reading grid inner intersections
    grid_intersections = gpd.read_file(grid_intersections_fp)
    
    # print(grid_intersections.head())
    
    # Preparing input to Pool:
        # the subset of inner intersections labelled with the same cell_id
        # subset of files with matching code: containing one of the grid_ids for the corresponding cellID
        # the identificator to write the files with the merged polygons (eg cID_16<cell16ID>_cID_32<cell32ID>)
    # poolInput = [(grid_intersections.loc[grid_intersections[grid_id] == cID, [grid_id, 'geometry']],
    #               fileList.selecting_fp_withCondition(files, [previous_grid_id, grid_id + str(cID)]),
    #               grid_id + str(cID) + '_' + next_grid_id + str(initial_grid.loc[initial_grid[grid_id] == cID, next_grid_id].tolist()[0]),
    #               year) for cID in cell_id]
    # poolInput = [(grid_intersections.loc[grid_intersections[next_grid_id] == cID, [next_grid_id, 'geometry']],
    #               fileList.selecting_fp_withCondition(files, [grid_id + str(e) for e in initial_grid.loc[initial_grid[next_grid_id] == cID, grid_id]])
                  
    #     ) for cID in cell_id]
    
    # print('Creating pool input')
    # Groups of inner intersections based on their next_grid_id id
    groups_grid_intersections = [grid_intersections[grid_intersections[next_grid_id] == cID] for cID in cell_id]
    
    # Groups of files containing polygons belonging to the same cell in next_grid_id
    groups_polygon_files = []
    for cID in cell_id:
        # Subset those cells in grid_id that correspond to cell_id in next_grid_id
        subset_cells_id = initial_grid.loc[initial_grid[next_grid_id] == cID, grid_id].tolist()
        # IMPORTANT: there was an important error here. I was missing requiring '_' after cID_<><> and so
        #            the algorithm was selecting a bunch of files for merging that shouldn't have been selected.
        #            For example, for cID_0820, it was also picking up cID_08200, cID_08204, etc.
        subset_cells_id = [grid_id + str(element) + '_' for element in subset_cells_id]
        
        subset_files_cells_id = fileList.selecting_fp_withCondition(files, subset_cells_id)
        # print(cID)
        # print(subset_cells_id)
        # print(subset_files_cells_id)
        
        # Append list of files 
        groups_polygon_files.append(subset_files_cells_id)
    # print('List of files to read form and then erase')   
    # print(groups_polygon_files)
    # Create list of codes for the filenames of the next grid
    next_grid_filecode = [next_grid_id + str(cID) for cID in cell_id]
    # print(next_grid_filecode)
    
    poolInput = list(zip( groups_grid_intersections, groups_polygon_files, next_grid_filecode, 
                         [year]*len(groups_polygon_files),
                         [output_polygon_folder]*len(groups_polygon_files)))
    
    
    with Pool(num_processes) as p:
        print('Starting parallel computation within merging_polygons_in_grid.')
        p.map(merge_polygons, poolInput)
        print('Finished parallel computation within merging_polygons_in_grid.')
        
        
    
    
def merge_polygons_all(initial_grid, grids_info_fp, year, num_process):
    """
    Function that gets the initial grid, the information on which grids to use and
    in what order, along with the folder from which to read and write to for each
    grid, along with the year in which the program is working and the number
    of processes to use in the parallelization.

    Parameters
    ----------
    initial_grid : GeoDataFrame
        The initial grid in which the polyon data is distributed.
    grids_info_fp : string
        Filepath to the CSV file containing the information on which to grids to
        use, their order and the input and output files for each one.
    year : int
        Year of data we are working on.
    num_process : int
        Number of processes to use in parallelization.

    Returns
    -------
    None.

    """
    
    # Reading the information on the different grids to use and their inner intersections
    grids_info = pd.read_csv(grids_info_fp)
    
    # We run the merging of polygons consecutively by grid level starting from the 
    # smaller sized cell grid. We do not do this for the last grid (a square) as 
    # there won't be any more polygons to merge
    for i in range(len(grids_info) - 1):
        
        t.tic()
        
        # Reading the grid ID
        grid_id = grids_info.loc[i, 'grid_id']
        
        print('\n Working at level grid {}'.format(grid_id))
        
        # Folder from whic to read the polygons
        input_polygon_folder = grids_info.loc[i, 'input_polygons_folder']
        # Folder where to write the output polygons
        output_polygon_folder = grids_info.loc[i, 'output_polygons_folder']
        
        # Reading filepath to the crosses of this grid
        grid_intersections_fp = grids_info.loc[i, 'grid_intersections']
        
        # Reading next grid ID
        next_grid_id = grids_info.loc[i + 1, 'grid_id']
        
        # Calling  function to merge the polygons that are touching one antoher 
        # in current grid
        merging_polygons_in_grid(initial_grid, input_polygon_folder, output_polygon_folder, grid_intersections_fp, grid_id, next_grid_id, year, num_process)

        # # If we are not working on the first grid, erase the temporary files
        # if i > 0:

        # Erase temporary files
        # First, get list of the files in folder temp
        files_to_erase = fileList.fileList(input_polygon_folder.format(year))
        # Subset those containing the grid_id
        files_to_erase = fileList.selecting_fp_withCondition(files_to_erase, [grid_id])

        for filepath in files_to_erase:
            os.remove(filepath)
        
        t.toc('Merged files at level {} in'.format(grid_id))
        
        

"""
MAIN
"""
if __name__ == "__main__":
    
    #---------------------- User inputs -------------------------------
    # Asking the user whether to read server orders for a specific node
    server_order = input("Input '_node<>' to select a specific set of orders.\nFor no-selection, just press enter.\n")
    
    # CSV with the list of raster files to work with in this script
    csv_fp = '../server_orders/MBfogo_raster_to_polygon{}.csv'.format(server_order)   
    
    # Output directory
    # out_dir = '../data/processed/MBfogo_c10_cerrado_polygons_{}.shp'
    
    # List of raster values to keep as polygons
    raster_values = list(range(1, 12 + 1))
    
    # Type of connectivity to be used
    connectivity = 8
    
    # Number fo processes for parallelization
    num_process = 38
    
    # #---------------------- Testing -----------------------------------
    # Read grid file of 0.8 degrees covering the extent of MapBiomas fogo
    # collection 1.0 over the Cerrado
    grid_fp = '../../data/shapes_ch/cerrado_grid_02deg_ch_cropped_labelled.shp'
    # grid_fp = '../data/overlays/cerrado_grid_08deg_cropped_labelled.shp'
    # grid_fp = '../data/overlays/cerrado_grid_08deg_top-right-quadrant.shp'
    # grid_fp = '../data/overlays/cerrado_grid_08deg_bottom-left-quadrant.shp'
    grid = gpd.read_file(grid_fp)
    
    #---------------------- The programme -----------------------------
    
    # Reading list of raster filepaths
    list_raster_fp = fileList.readListFiles_fromCSV(csv_fp)
      
    
    
    # Polygonise each raster file in the list
    for in_fp in list_raster_fp:
        print("Working on file {}".format(os.path.basename(in_fp)))
        
        # Extract year of data from filepath
        year = extract_time.extract_year(in_fp)
        
        # # ----------------- Merging polygons ----------------------------------
        # # Merge polygons by window and month
        
        # # Drop geometry in grid table as no longer needed
        # grid = grid.drop(columns = ['geometry'])
        
        
        print('Starting merging polygons by grid level.')
        merge_polygons_all(grid, '../aid_inputs/raster_to_polygon_merging_info_c20_1999.csv', int(year), num_process)


        
    

