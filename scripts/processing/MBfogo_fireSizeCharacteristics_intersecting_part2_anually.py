# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 11:10:26 2023


@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

**IMPORTANT: This script must be run in the scipy-env environment.

Script that calculates the fire size characteristics per cell and period. The 
fire characteristics calculated are:
    - 5% quantile
    - 50% quantile (median)
    - 95% quantile
    - the exponent (alpha) and sigma (standard error) of a fitted power law
    
    
Input files:
    - ../../data/processed/MBfogo_c{}0_cerrado_polygons_b150_CM_fixed/MBfogo_c{}0_cerrado_polygons_{}.csv

Output files:
    - ../../data/processed/summary_tables/fireSize_intersecting_grid{}_period{}.csv'  
"""

"""
Imports
"""
import pandas as pd
import geopandas as gpd
from multiprocessing import Pool
from pytictoc import TicToc
import numpy as np
t = TicToc()

import powerlaw

import sys
sys.path.append('../utils/')
import parallel_prep as pprep


"""
Functions
"""
# def collectArea_intersectingPolygon(slices, grid_fp, polygon_fp, polygon_data_fp, min_size):

#     # Read chunk of cells in grid
#     cells = gpd.read_file(grid_fp, rows = slice(slices[0], slices[1]))
    
#     # Convert cells chunk to poylgons CRS
#     cells = cells.to_crs('epsg:4326')
    
#     # # Read chunk of polygon data in csv file
#     # polygons_data = pd.read_csv(polygon_data_fp, names = col_names, skiprows = lambda x: x not in range(slices[0] + 1, slices[1]+1+1))
#     # # print(polygons_data)
#     # # Convert grid CRS to that of polygon
#     # grid = grid.to_crs(polygons.crs)
    
#     # Creating empty lists to accumulate the areas of the polygons that intersect the cells
#     df_area_cell= []
    
#     # Loop over cells to find the intersecting polygons in year
#     for i in range(len(cells)): 
        
#         # Cell ID
#         cell_id = cells.loc[i, 'id'].item()
        
#         # Read intersecting poylgons
#         int_polygons = gpd.read_file(polygon_fp, bbox = cells.loc[i, 'geometry'])
        
#         # Selecting only the polygon id
#         pol_id = int_polygons['ID'].tolist()
        
#         # Reading the polygons data
#         pol_data = pd.read_csv(polygon_data_fp)
        
#         # Selecting the areas of the intersecting polygons that are larger than 0.03 km2
#         pol_areas = pol_data.loc[(pol_data['ID'].isin(pol_id)) & (pol_data['area'] >= min_size), 'area'].tolist()
        
#         # Creating dataframe with the information
#         # id (cell) | area
#         cell_areas = pd.DataFrame({'id' : cell_id, 'area': pol_areas})
        
#         # Append to list
#         df_area_cell.append(cell_areas)
        
#     # Cast list of data frames as data frame
#     df_area_cell = pd.concat(df_area_cell, ignore_index = True)
    
#     return df_area_cell



def getQuantiles(df, list_q):
    """
    Function that calculates the quantiles of fire size (area) indicated in the 
    list of quantiles, list_q, for each cell in dataframe. Returns a dataframe
    of the form:
        id | q_<>

    Parameters
    ----------
    df : DataFrame
        Dataframe containing a list of fire  polygons, each one with the grid cell
        they belong to (id) and their area in km2 (area).
    list_q : List of floats
        The different quantiles we want to describe the fire size distribution in each cell.

    Returns
    -------
    df_cell_quantiles : DataFrame
        Dataframe of the form: id | q_<>

    """
    
    # Grouping by cell and getting the different quantiles
    df_cell_quantiles = df[['id', 'area']].groupby(by = 'id').quantile(q = list_q).reset_index()
    
    # Casting from wide to long
    df_cell_quantiles = df_cell_quantiles.pivot(index = 'id', columns = 'level_1', values = 'area').reset_index()
    
    df_cell_quantiles.columns = ['id'] + ['q_'+ str(q) + '_int' for q in list_q]
    
    return df_cell_quantiles


def fitPowerlaw(df, x_min):
    """
    Function that fits a powerlaw to a fire size distribution. Returns teh input
    dataframe with two additional columns indicating the the exponent of the
    distribution and its standard error.

    Parameters
    ----------
    df : DataFrame
        Dataframe where each row is a different polygon and contains the columns
        'id' indicating the cell the polygon belongs to, and its 'area' in km2.
    x_min : float
        The minimu possible size of the fire polygon size distribution.

    Returns
    -------
    df : DataFrame
        Dataframe where each row is a cell, with the fitted power-law exponent
        alpha and its standard error. Of the form: id | pwl_alpha | pwl_sigma.

    """

    # Group by cell and create powerlaw.Fit object for each group of areas, 
    # passing the xmin value provided by user
    df = df[['id', 'area']].groupby(by = 'id').agg({'area': lambda x: powerlaw.Fit(x, xmin = x_min)})

    
    # Changing column name for clarity
    df = df.rename(columns = {'area': 'fit'})
    
    # Calculating the exponent (alpha) and its standard error (sigma)
    df['pwl_alpha_int'] = df['fit'].apply(lambda x: x.power_law.alpha)
    df['pwl_sigma_int'] = df['fit'].apply(lambda x: x.power_law.sigma)
    df = df.drop(columns = 'fit').reset_index()

    
    return df

"""
MAIN
"""
def main():
    
    # ---------------- User inputs --------------------------------------------
    
    # Grid to work on 
    grid_id = '30km'
    grid_col_name = 'cellID_{}'.format(grid_id)
    
    
    # Start and end years
    start_year = 1985
    end_year = 2020
    
    # Choose the MapBiomas FOGO collection
    collection = 2
    
    # Filepath to files
    cell_size_period_fp = '../../data_new_grids/processed/summary_tables/annual/fireSize_intersecting_perCell_grid{}_year{}.csv'
    
    # Output file
    output_fp = '../../data_new_grids/processed/summary_tables/annual/fireSize_intersecting_grid{}_year{}.csv'
    
    
    # Minimum fire size
    min_polSize = 0.03

    # list of quantiles we want the information for
    list_quantiles = [0.05, 0.5, 0.95, 0.99, 1.]
    
    #--------------------------------------------------------------------------
    
    # -------------- Fire characteristics -------------------------------------
    
    # We work on a group of years at a time
    for year in range(start_year, end_year + 1):
        
        print('Working on year {}'.format(year))
        
        # Reading the data, which comes in the format: id (cell) | period | area
        data_year = pd.read_csv(cell_size_period_fp.format(grid_id, year))
        
        # Calculating area quantiles for each cell
        data_quantiles = getQuantiles(data_year, list_quantiles)
        
        # Fitting a powerlaw to each cell's areas and calculating the exponent
        # the corresponding standard error
        data_pwl = fitPowerlaw(data_year, min_polSize)
        
        # Merging the two bits of data
        data_year = data_quantiles.merge(data_pwl, on = 'id', how = 'outer')
        
        # Freeing memory
        del data_quantiles, data_pwl
        
        # Adding column indicating the period
        data_year['year'] = year
        
        
        # Reordering columns
        data_year = data_year[['id', 'year'] + data_year.columns[1:-1].tolist()]
        
        # print(data_period)
        
        # Saving to file
        data_year.to_csv(output_fp.format(grid_id, year), index = False)
        
        
        print('')
    
    
    
    
if __name__ == '__main__':
    main()