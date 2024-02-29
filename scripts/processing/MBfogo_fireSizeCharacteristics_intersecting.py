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
t = TicToc()

# import powerlaw

import sys
sys.path.append('../utils/')
import parallel_prep as pprep


"""
Functions
"""
def collectArea_intersectingPolygon(slices, grid_fp, polygon_fp, polygon_data_fp, min_size):

    # Read chunk of cells in grid
    cells = gpd.read_file(grid_fp, rows = slice(slices[0], slices[1]))
    
    # Convert cells chunk to poylgons CRS
    cells = cells.to_crs('epsg:4326')
    
    # # Read chunk of polygon data in csv file
    # polygons_data = pd.read_csv(polygon_data_fp, names = col_names, skiprows = lambda x: x not in range(slices[0] + 1, slices[1]+1+1))
    # # print(polygons_data)
    # # Convert grid CRS to that of polygon
    # grid = grid.to_crs(polygons.crs)
    
    # Creating empty lists to accumulate the areas of the polygons that intersect the cells
    df_area_cell= []
    
    # Loop over cells to find the intersecting polygons in year
    for i in range(len(cells)): 
        
        # Cell ID
        cell_id = cells.loc[i, 'id'].item()
        
        # Read intersecting poylgons
        int_polygons = gpd.read_file(polygon_fp, bbox = cells.loc[i, 'geometry'])
        
        # Selecting only the polygon id
        pol_id = int_polygons['ID'].tolist()
        
        # Reading the polygons data
        pol_data = pd.read_csv(polygon_data_fp)
        
        # Selecting the areas of the intersecting polygons that are larger than 0.03 km2
        pol_areas = pol_data.loc[(pol_data['ID'].isin(pol_id)) & (pol_data['area'] >= min_size), 'area'].tolist()
        
        # Creating dataframe with the information
        # id (cell) | area
        cell_areas = pd.DataFrame({'id' : cell_id, 'area': pol_areas})
        
        # Append to list
        df_area_cell.append(cell_areas)
        
    # Cast list of data frames as data frame
    df_area_cell = pd.concat(df_area_cell, ignore_index = True)
    
    return df_area_cell



# def getQuantiles(df, list_q):
#     """
#     Function that calculates the quantiles of fire size (area) indicated in the 
#     list of quantiles, list_q, for each cell in dataframe. Returns a dataframe
#     of the form:
#         id | q_<>

#     Parameters
#     ----------
#     df : DataFrame
#         Dataframe containing a list of fire  polygons, each one with the grid cell
#         they belong to (id) and their area in km2 (area).
#     list_q : List of floats
#         The different quantiles we want to describe the fire size distribution in each cell.

#     Returns
#     -------
#     df_cell_quantiles : DataFrame
#         Dataframe of the form: id | q_<>

#     """
    
#     # Grouping by cell and getting the different quantiles
#     df_cell_quantiles = df[['id', 'area']].groupby(by = 'id').quantile(q = list_q).reset_index()
    
#     # Casting from wide to long
#     df_cell_quantiles = df_cell_quantiles.pivot(index = 'id', columns = 'level_1', values = 'area').reset_index()
    
#     df_cell_quantiles.columns = ['id'] + ['q_'+ str(q) for q in list_q]
    
#     return df_cell_quantiles


# def fitPowerlaw(df, x_min):
#     """
#     Function that fits a powerlaw to a fire size distribution. Returns teh input
#     dataframe with two additional columns indicating the the exponent of the
#     distribution and its standard error.

#     Parameters
#     ----------
#     df : DataFrame
#         Dataframe where each row is a different polygon and contains the columns
#         'id' indicating the cell the polygon belongs to, and its 'area' in km2.
#     x_min : float
#         The minimu possible size of the fire polygon size distribution.

#     Returns
#     -------
#     df : DataFrame
#         Dataframe where each row is a cell, with the fitted power-law exponent
#         alpha and its standard error. Of the form: id | pwl_alpha | pwl_sigma.

#     """
    
#     # Group by cell and create powerlaw.Fit object for each group of areas, 
#     # passing the xmin value provided by user
#     df = df[['id', 'area']].groupby(by = 'id').agg({'area': lambda x: powerlaw.Fit(x, xmin = x_min)})
    
#     # Changing column name for clarity
#     df = df.rename(columns = {'area': 'fit'})
    
#     # Calculating the exponent (alpha) and its standard error (sigma)
#     df['pwl_alpha'] = df['fit'].apply(lambda x: x.power_law.alpha)
#     df['pwl_sigma'] = df['fit'].apply(lambda x: x.power_law.sigma)
#     df = df.drop(columns = 'fit').reset_index()
    
#     return df

"""
MAIN
"""
def main():
    
    # ---------------- User inputs --------------------------------------------
    
    # Grid to work on 
    grid_id = '50km'
    grid_col_name = 'cellID_{}'.format(grid_id)
    
    # Number of periods
    n_periods = 4
    
    if n_periods == 2:
        # Year periods to work on 
        periods = [(1985, 2002), (2003, 2020)]
        
        # Number of years per period
        n_years = 18
        
    elif n_periods == 4:
        # Year periods to work on 
        periods = [(1985,1993), (1994, 2002), (2003, 2011), (2012, 2020)]
        
        # Number of years per period
        n_years = 9
    
    # period_labels = [str(p) for p in periods]
    
    # # Start and end years
    # start_year = 1985
    # end_year = 2020
    
    # UPDATE (24/10/2023): Allowing to choose the MapBiomas FOGO collection
    collection = 2
    
    # Filepath to files
    data_year_fp = '../../data/processed/MBfogo_c{}0_cerrado_polygons_b150_CM/MBfogo_c{}0_cerrado_polygons_{}.{}'
    
    # Output file
    output_fp = '../../data/processed/summary_tables/{}_year_periods/fireSize_assigned_perCell_grid{}_period{}.csv'


    # UPDATE (24/10/2023): We calculate the data for all grid cells
    # Work only with cell that have more than a certain percentage of their area inside the Cerrado
    # thresh = 0.0 # 75.0
    
    # # Extra cells to include and cells to exclude
    # include_cells = [567, 568]
    # exclude_cells = [1451]
    
    # List of quantiles we want to work on
    # list_quantiles = [0.05, 0.5, 0.95, 0.99, 1.]
    
    # The minimum polygon area that we accept
    min_polSize = 0.03
    
    num_process = 38
    
    
    #--------------------------------------------------------------------------
    
    # ------------- Reading data ----------------------------------------------
    # print('Reading data for all years...\n')
    
    # # Read and format the data for all years consecutively
    # data_year = readFormat_dataYears(data_year_fp, start_year, end_year, grid_id, ['area'], 
    #                                  periods = True, nperiods = n_periods, grid_colname = grid_col_name, 
    #                                  fire_polygons = True, min_size = min_polSize, 
    #                                  thresh = thresh)
    
    # print(data_year)#
    # print(data_year.columns)
    
    # Reding grid
    grid_fp = '../../data/shapes/cerrado_grid_{}_cropped.shp'.format(grid_id)
    grid = gpd.read_file('../../data/shapes/cerrado_grid_{}_cropped.shp'.format(grid_id))
    # # grid = grid[['id','geometries']]
    # print(grid)
    # print(len(grid))
    # -------------------------------------------------------------------------
    
    # -------------- Fire characteristics -------------------------------------
    
    i = 1
    # We work on a group of years at a time
    for (start_year, end_year) in periods:
        
        period = '({}, {})'.format(start_year, end_year)
        
        print('Working on group ' + period)
        
        # Creating empty lists where we accumulate the sizes of all intersecting polygons per cell
        df_sizes = []
        
        # For each year in the period, 
        for year in range(start_year, end_year + 1):
            
            print('')
            print('Working on year {}'.format(year))
        
            # ---------------------------------------------------------------------
            # Getting paths to input files
        
            # Filepath to files
            files = data_year_fp.format(collection, collection, '{}', '{}')
            # Filepath to fire polygon shapefile
            polygon_fp =  files.format(year, 'shp')
            
            # Filepath to fire polygon data
            polygon_data_fp = files.format(year, 'csv')
            
            # Parallelise over the cells so find all intersecting polygons
            
            # Only the total number of rows is needed to be able to produce partitions
            ncells = len(grid)
            
            # Generate the slices
            slices = pprep.list_slices(ncells, num_process)
            # print(slices)
            
            # Generating input to Pool using the slices
            pool_input = list(zip(slices, 
                                  [grid_fp]*len(slices),
                                  [polygon_fp]*len(slices), 
                                  [polygon_data_fp]*len(slices),
                                  [min_polSize]*len(slices)))
            
            with Pool(num_process) as p: 
                
                print("Starting parallel computation in collectArea_intersectingPolygon...")
                
                t.tic()
                df_cell_area = p.starmap(collectArea_intersectingPolygon, pool_input)

                t.toc('Finished parallel computation.')
                
                
            # Converting list of dataframes to a single dataframe
            df_cell_area = pd.concat(df_cell_area, ignore_index = True)
            # print(df_cell_area.head())
            
            # Appending to list of dataframes with sizes, one per year
            df_sizes.append(df_cell_area)
            
        # Cast list of dataframes as unique dataframe
        df_sizes = pd.concat(df_sizes, ignore_index = True)
        
        # Now we add a column to each indicating the period
        df_sizes['period'] = str(periods[i-1])
        
        # Reshuffling columns
        df_sizes = df_sizes[['id','period','area']]
        
        # Saving to file
        df_sizes.to_csv(output_fp.format(n_years, grid_id, i), index = False)
        i += 1
            
        
        
        # Calculating area quantiles for each cell
        # data_quantiles = getQuantiles(df_sizes, list_quantiles)
        
        # # Fitting a powerlaw to each cell's areas and calculating the exponent
        # # the corresponding standard error
        # data_pwl = fitPowerlaw(df_sizes, min_polSize)
        
        # # Merging the two bits of data
        # data_period = data_quantiles.merge(data_pwl, on = 'id', how = 'outer')
        
        # # Freeing memory
        # del data_quantiles, data_pwl
        
        # # Adding column indicating the period
        # data_period['period'] = periods[i]
        # i += 1
        
        # # Reordering columns
        # data_period = data_period[['id', 'period'] + data_period.columns[1:-1].tolist()]
        
        # print(data_period)
        
        # # Saving to file
        # data_period.to_csv(output_fp.format(n_years, grid_id, i), index = False)
        
        print('')
    
    
    
    
if __name__ == '__main__':
    main()