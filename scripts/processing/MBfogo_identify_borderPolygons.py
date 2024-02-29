# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 14:22:09 2022

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that identifies those polygons falling within 30 m (one Landsat pixel) 
of the border of the Cerrado.

Input files:
    ../../data/processed/MBfogo_c10_cerrado_polygons_b150_CM/MBfogo_c10_cerrado_polygons_<year>.<>
    ../../data/shapes/cerrado_1-250000_CRS5880_borderBuffer30.shp
    ../../data/shapes/cerrado_grid_<grid_size_units>_cropped.shp

Output files (overwriting input files):
    ../../data/processed/MBfogo_c10_cerrado_polygons_b150_CM/MBfogo_c10_cerrado_polygons_<year>.csv
    ID | year | month | ... | border_pol



This script has the following functions:
    * id_borderPolygons - finds the set of polygons that intersect the border
    * id_borderPolygons_withPool - prepares the data to run id_borderPolygons 
                                   in parallel
    * id_borderCandidates - identify polygons that are candidates to be at the 
                            border of the Cerrado




**Method to create the border polygon
border = gpd.read_file(border_fp)
# Extracting the boundary
border = border.boundary.buffer(30)#.plot(figsize=(20,15))
border.to_file('../../data/shapes/cerrado_1-250000_CRS5880_borderBuffer30.shp')

"""

"""
Imports
"""
import geopandas as gpd
import pandas as pd
from multiprocessing import Pool
from pytictoc import TicToc
t = TicToc()

import sys
sys.path.append('../utils/')

import parallel_prep as pprep
import fileList

"""
Functions
"""
def id_borderPolygons(geom_candidates, frontier):
    """
    Function that finds the set of polygons in geom_candidates that intersect 
    with the frontier

    Parameters
    ----------
    geom_candidates : GeoDataFrame
        GeoDataFrame containing the polygons that may intersect with the frontier.
    frontier : shapely geometry (polygon)
        The geometry of the border polygon.

    Returns
    -------
    geom_border : list of strings
        List of IDs of those poylons that intersect with frontier.

    """
    
    # Subset IDs of geometries intersecting the border
    geom_border = geom_candidates.loc[geom_candidates['geometry'].intersects(frontier), 'ID']
    
    # Return list of border geometries
    return geom_border
    



def id_borderPolygons_withPool(candidates, border, n_process):
    """
    Function that identifies the polygons that intersect with the border out 
    of all the candidates. It first prepares the data so that the function 
    id_borderPolygons can be run in parallel.

    Parameters
    ----------
    candidates : List of strings
        list of polygon IDs of all polygons that are candidates to being borer polygons.
    border : shapely geometry (polygon)
        The geometry of the border polygon.
    n_process : integer
        Number of processes to use in parallelisation.

    Returns
    -------
    list_polygon_border : list of strings
        List of polygon identifiers that are in the border (intersect with border).

    """

    print('Preparing Pool input. \n')
    
    # First divide this geodataframe of candidates into differnet chunks
    # to run id_borderPolygons in parallel
    slices = pprep.list_slices(len(candidates), n_process)
    
    # Use these slices to convert the dataframe into a list of dataframes
    # containing each one a slice of the original one
    list_candidates = []
    for chunk in slices:
        
        # For how slices work in comparison to Python index subsetting, we
        # have to add one unit to the last index so that we read all rows
        correction = 0
        if chunk == slices[-1]:
            correction = 1
        
        # Subset and append to the list
        list_candidates.append(candidates.iloc[chunk[0]:chunk[1] + correction, :])

    
    # Generating input to Pool using the slices
    pool_input = list(zip(list_candidates, [border]*len(slices)))
    

    with Pool(n_process) as p: 
        
        print("Starting parallel computation in id_borderPolygons_withPool...")
        
        list_polygon_border = p.starmap(id_borderPolygons, pool_input)
        
        print('Finished parallel computation. \n')

    
    # List of border polygons is a list of lists, casting as a list
    list_polygon_border = [item for sublist in list_polygon_border for item in sublist]
    
    return list_polygon_border
    



def id_borderCandidates(polygon_csv_fp, polygon_shp_fp, grid_size_units, list_border_cells):
    """
    Function that identifies border candidate polygons: those polygons found within
    a bordering grid cell. Returnsthe list of IDs of the candidate polygons

    Parameters
    ----------
    polygon_csv_fp : string
        Filepath to CSV polygon dataset with polygon IDs and other information. 
        One of the columns must be called callID_<grid_size_units> and indicate 
        the cell ID to which the polygon belongs.
    polygon_shp_fp : string
        Filepath to SHP polygon dataset containing the polygons to be identified.
    grid_size_units : string
        Grid identifier with the grid size and the units used. E.g. 50km.
    list_border_cells : list of integers
        List containing the identifiers of the grid cell that are found in the 
        border.

    Returns
    -------
    polygon_candidates : list of string
        List of identifiers of the polygons that are candidates to being border 
        polygons.

    """
    
    # First, read csv datafile where to add column indicating if border polygon
    polygon_data = pd.read_csv(polygon_csv_fp)
    
    
    # Using grid to subset candidates to border polygons - getting just the ID
    border_candidates = polygon_data.loc[polygon_data['cellID_{}'.format(grid_size_units)].isin(list_border_cells), 'ID']


    # Freeing memory - we don't need this dataframe for the moment
    del polygon_data
    
    
    print('Reading polygons to identify cnadidates...')
    # Now we read the shapefile containing all the poylgons, subset the ones 
    # we are interested and erase the rest
    polygons = gpd.read_file(polygon_shp_fp)
    

    print('Subsetting polygon candidates. \n' )
    polygon_candidates = polygons[polygons['ID'].isin(border_candidates)].reset_index(drop = True)
    
    # Freeing memory 
    del polygons
    
    return polygon_candidates
    
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
    
    
    # CSV with the list of raster files to work with in this script
    csv_fp = '../server_orders/MBfogo_identify_borderPolygons{}.csv'.format(server_order) 
    
    # Reading list of files to work on (one after the other)
    list_years = fileList.readListFiles_fromCSV(csv_fp)

    # Grid size
    grid_size_units = list_years[0]
    
    
    # The number of processes to parallelise over
    num_process = 35
    
    
    # Input filepaths
    
    # Generic filepath to polygon files (input files)
    filepath = '../../data/processed/MBfogo_c10_cerrado_polygons_b150_CM/MBfogo_c10_cerrado_polygons_{}.{}'
    
    # File containing polygon whose border we want to check against
    border_fp = '../../data/shapes/cerrado_1-250000_CRS5880_borderBuffer30.shp'
    
    # Filepath to grid - we need the list of border cells
    grid_fp = '../../data/shapes/cerrado_grid_{}_cropped.shp'.format(grid_size_units)
    
    

    # -------------- Identifying border polygons by year ----------------------
    
        
    # Reading border geometry (already with a 30m buffer on both sides)
    border = gpd.read_file(border_fp)
    
    
    # Subset list of cells not fully within the Cerrado
    # Reading grid
    grid = gpd.read_file(grid_fp)
    # Getting list of cells at the border of the Cerrado (not fully within the Cerrado)
    border_cells = grid.loc[grid['withinCerr'] == 0, 'id']
    
    
    # Freeing memory
    del grid
    
    # Converting list of string to list of integers
    list_years = [int(year) for year in list_years[1:]]
    
    
    # For each year, identify and label border polygons
    for year in list_years:
        
        t.tic()
        print('')
        print('Working on year {}'.format(year))
        
        # ------------------- Datafiles to use -------------------------------
        
        # Filepath to fire polygon shapefile
        polygons_fp = filepath.format(year, 'shp')
        
        # Filepath to fire polygon data
        polygon_data_fp = filepath.format(year, 'csv')
        
        # --------------------------------------------------------------------
        
        
        # ------------- Labelling border polygons ----------------------------
        
        # 1. Identify polygon candidates to being border polygons
        #    i.e. those that fall in border cells
        polygon_candidates = id_borderCandidates(polygon_data_fp, polygons_fp, grid_size_units, border_cells)
        

        # Convert the border polygon CRS to that of the fire polygon data and 
        # extract the geometry. Only needed once.
        if year == list_years[0]:
            # Convert CRS
            border = border.to_crs(polygon_candidates.crs)
            # and then we extract its geometry
            border = border.loc[0, 'geometry']
        
        print('\nFinding candidates to border polygons. \n')
        

        # 2. Identify border polygons out of candidates 
        #    i.e.those that intersect with the border geometry
        list_polygon_border = id_borderPolygons_withPool(polygon_candidates, border, num_process)
        
        
        # Freeing memory
        del polygon_candidates
        
        
        # 3. Finally, record it in CSV
        # Using the list of IDs of border polygons
        # Read csv datafile again
        print('Reading polygon data and labelling border polygons')
        polygon_data = pd.read_csv(polygon_data_fp)
        
        # Create column to record border polygons. Initialise to 0 (not border polygon)
        polygon_data['border_pol'] = 0
        
        # Labelling border polygons
        polygon_data.loc[polygon_data['ID'].isin(list_polygon_border), 'border_pol'] = 1
        
        
        # ------------- Saving output ----------------------------------------
        
        # Writing file again with this information
        print('\n Writing final dataframe')
        polygon_data.to_csv(polygon_data_fp, index = False)
        #---------------------------------------------------------------------
        
        t.toc('Year {} done in'.format(year))
        
        print('\n \n')


if __name__ == '__main__':
    main()
