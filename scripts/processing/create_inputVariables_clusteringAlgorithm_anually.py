# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 10:19:31 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

UPDATE 24/08/2023: Change in the variables it considers to include only the 
                   ones to be used in the clustering algorithm. As well, updated
                   input and output files to new system in which different sets
                   of periods are allowed.

Script that creates a unique summary of variables per period. Variables are the
input candidates for the clustering algorithm.

Variables selected:
    - median of the fire size
    - 99th of the fire size
    - exponent of the power law
    - number of assigned fires
    - 99th of the frequency distribution
    - median of the season (excluding the 3 rainy season months and adjusting 
                            for the fire year)
    - burned area (of the polygons intersection)
    
Input files:
    ../../data/processed/summary_tables/{}_year_periods/
    - fireSize_assigned_grid{}_period{}.csv
    - burnedArea_fireCount_assigned_grid{}_period{}.csv
    - fireFrequency_returnInterval_grid{}_period{}.csv
    - fireSeason_measures_exclude3RainyMonths_grid{}_period{}.csv
    - burnedArea_fireCount_intersection_monthFromPolygon_grid{}_period{}.csv
    
Output file:
    ../../data/processed/clustering_algorithm/{}_year_periods/input_vars_grid{}_period{}.csv
    

"""

"""
Imports
"""
import pandas as pd
import geopandas as gpd
import numpy as np
from functools import reduce


"""
Functions
"""
def replace(my_list, my_dict):
    return [x if x not in my_dict else my_dict[x] for x in my_list]


def spatial_average(cell, data, var, edge_list):
    # Function that, given a cell and a set of data, and a variable, and a connectivity list, 
    # calculates the average
    
    # Fetching the list of cells connected to cell
    list_cells = edge_list[(edge_list['from'] == cell) | (edge_list['to'] == cell)][['from', 'to']]
    list_cells = list_cells['from'].tolist() + list_cells['to'].tolist()
    list_cells = [element for element in list_cells if element != cell]

    
    # Calculating the average value
    spatial_average = data.loc[data['id'].isin(list_cells), var].mean()
    
    
    return spatial_average
    


"""
MAIN
"""
def main():
    
    # -------------------------- User inputs ----------------------------------
    
    # Grid id
    grid = '30km'
    
    # Study period
    start_year = 1985
    end_year = 2020
    
    # Files
    #--------
    
    # Input directory
    in_dir = '../../data_new_grids/processed/summary_tables/annual/'
    
    # Filename end
    file_termination = '_grid{}_year{}.csv'
    
    # File identifiers
    # fire_size_fp = 'fireSize_assigned'                                        # UPDATE 31/10/2023 Method dismissed in favour of intersecting
    number_fires_fp = 'burnedArea_fireCount_assigned'                           # All cells
    frequency_fp = 'fireFrequency_fromPixels_sECDF_quantiles'                   # UPDATE 26/10/2023 Including it again as it may be less correlated with burned area
                                                                                # UPDATE 24/08/2023  # 'fireFrequency_fromPixels_quantiles'        
    fire_return_interval_fp = 'fireFrequency_returnInterval'                 
    season_fp = 'fireSeason_measures_exclude3RainyMonths'
    burned_area_fp = 'burnedArea_fireCount_intersection_monthFromPolygon'
    fire_size_fp = 'fireSize_intersecting'                                  # UPDATE 31/10/2023 Different calculation of the fire size by taking all intersecting fires instead of the assigned ones
    
    # Cast file identifiers into list
    file_ids = [fire_size_fp, number_fires_fp, season_fp, burned_area_fp]
    
    
    # Output filepath
    out_fp = '../../data_new_grids/processed/clustering_algorithm/annual/input_vars_grid{}_year{}.csv'
    
    
    # Variables
    #-----------
    fire_season_var = ['median_month_sECDF', 'median_month_sECDF_rel']          # ['ci', 'median_month'] UPDATE 24/08/2023
    # fire_size = ['q_0.5', 'q_0.99', 'pwl_alpha']                              # UPDATE 31/10/2023 Method dismissed in favour of intersecting
    number_fires_var = ['npolygons']
    burned_area_var = ['area_T']
    # frequency_var = ['f99']                                                     # UPDATE 26/10/2023 including UPDATE 24/08/2023
    # fri_var = ['fri']
    # fri_freq_var = ['fri_freq']                                                 # UPDATE 26/10/2023 including 
    flammable = ['flammable']                                                   # UPDATE(26/10/2023) including to normalise burned area and number of fires
    fire_size = ['q_0.5_int', 'q_0.99_int']                                # UPDATE 31/10/2023 Method intersecting favoured instead of assigned
    
    # Cast variables into list
    list_vars = fire_season_var + fire_size + number_fires_var \
        + burned_area_var + flammable                                           # + fri_var + fri_freq_var + frequency_var 

    list_vars_labels = ['season_q50s', 'season_q50s_rel', 'number_fires_n', 'burned_area_n',
                        'size_q50', 'size_q99']
    list_vars_transform = ['none', 'none', 'log', 'log', 'log', 'log']
    
    # Columns names
    #---------------
    
    # Adjusting column names to desired ones
    new_col_names = {'ci': 'season_ci',
                     # 'median_month': 'season_q50',                            UPDATE 24/08/2023
                     'median_month_sECDF': 'season_q50s',
                     'median_month_sECDF_rel' : 'season_q50s_rel',
                     'q_0.5_int': 'size_q50',
                     'q_0.99_int': 'size_q99',
                     # 'pwl_alpha': 'size_alpha',
                     'npolygons': 'number_fires',
                     # 'f99': 'frequency_q99s',
                     'area_T': 'burned_area'}                                  
    
    
    # UPDATE (26/10/2023)
    # Assigning mean neighbour value to seasonality measures of cells
    imputation_cells = {
        2376 : [2286, 2287, 2288, 2375, 2377, 2464, 2465, 2466], 
        2385 : [2295, 2296, 2297, 2384, 2386, 2473, 2474, 2475]}

    # grid
    #---------------
    grid_dat = gpd.read_file('../../data/shapes/cerrado_grid_{}_cropped.shp'.format(grid))
    # Keeping only columns of interest - %area_crop is only a dummy column
    grid_dat = grid_dat[['id','%area_crop']]
    
    
    # Edge list
    #---------------
    edge_list = pd.read_csv('../../data_new_grids/processed/clustering_algorithm/edge_list_connectivity8_grid_{}_epsg5880.csv'.format(grid))
    
    # ------------------ Reading and formatting input data --------------------
    
    # Working on one period at a time
    for year in range(start_year, end_year + 1):
        
        print('Working on year {}'.format(year))
        
        data = []
        
        # Reading files consecutively
        for file_fp in file_ids:
            print('Working on {}'.format(file_fp))
            data_file = pd.read_csv(in_dir + 
                                    file_fp + 
                                    file_termination.format(grid, year))
            
            # print(data_file)
            
            # Manual cleaning to avoid same-column names. 
            if file_fp == number_fires_fp:
                data_file = data_file.drop(columns = 'area_T')
                
            elif file_fp == burned_area_fp:
                data_file = data_file.drop(columns = 'npolygons')
            
            elif file_fp == fire_return_interval_fp:
                data_file = data_file.drop(columns = ['area_T'])                # UPDATE 26/10/2023 , 'flammable'])
            
            # Selecting relevant columns only
            data_file = data_file[data_file.columns.intersection(['id', 'year'] + list_vars)]

            
            # Adjusting names of columns
            data_file.columns = replace(data_file.columns, new_col_names)
            
            # print(data_file.head(10))
            
            # Appending to list of dataframes
            data.append(data_file)
            
            # Freeing memory
            del data_file
            
        # Reading land-use data to fetch flammable area
        data_file = pd.read_csv('../../data_new_grids/processed/lulc_data/MBlanduse_grid{}_cerrado_{}.csv'.format(grid, year))
        data_file['flammable'] = data_file[['1', '11', '12', '13', '14']].sum(axis = 1)
        data_file['year'] = year
        data_file = data_file[['id', 'year', 'flammable']]
        
        data.append(data_file)
        
        # print(data)
        # Casting list of dataframes to unique dataframe
        data = reduce(lambda x, y: x.merge(y, on = ['id','year']), data)
        
        
        # print(data)
        
        # UPDATE 26/10/2023
        # Right now we have all the data loaded for all the periods. But right here
        # there are cells that do not appear in the dataframe because they do not
        # contain any fires in this period.
        # We are going to perform data imputation for these cells so that they
        # have data and can be included in the PCA.
        # The data imputation is goign to work like this:
        # we assign one fire to the area with the smallest size possible (a fire 
        # of 3 ha or 0.03 km2), and the frequency is effectively 0
            # number_fires = 1
            # burned_area = 0.03
            # size_q50 = 0.03
            # size_q99 = 0.03
            # frequency_q99s = 0
            # flammable = spatial average
            # fri = spatial average
            # fri_freq = spatial average
            # season_q50s = spatial average
        
        # Similarly, for those cells that only have one fire (seasonality calculations
        # are not possible), we assign them the spatial average
        
        # First, we merge the grid data, which contains the full list of cells
        data = grid_dat.merge(data, on = 'id', how = 'left')
        # Dropping dummy variable
        data = data.drop(columns = '%area_crop')
        
        # UPDATE 
        
        
        # Fetching list of cells that do not have any fires
        empty_cells = data.loc[np.isnan(data['number_fires']), 'id'].tolist()
        
        if len(empty_cells) > 0:
            print('\nThere are empty cells...')
            
            # Filling period
            data.loc[data['id'].isin(empty_cells), 'year'] = year # data.loc[0, 'year']
        
            # Fill data accordingly
            data.loc[data['id'].isin(empty_cells), 'number_fires'] = 1
            data.loc[data['id'].isin(empty_cells), 'burned_area'] = 0.03
            data.loc[data['id'].isin(empty_cells), 'size_q50'] = 0.03
            data.loc[data['id'].isin(empty_cells), 'size_q99'] = 0.03
            # data.loc[data['id'].isin(empty_cells), 'frequency_q99s'] = 0
            data.loc[data['id'].isin(empty_cells), 'season_q50s'] = 0
            
            
            # Spatial averages
            for cell in empty_cells:
                print('Working on cell {}...'.format(cell))
                for var in ['flammable']:#, 'fri', 'fri_freq']: #,'size_alpha']:
                    data.loc[data['id'] == cell, var] = spatial_average(cell, data, var, edge_list)
        
                    
        
        # List of cells with only one fire
        oneFire_cells = data.loc[data['number_fires'] == 1, 'id'].tolist()

        while len(oneFire_cells) > 0:
            print('\nThere are cells with only one fire...')
        
            # Calculating the seasonality as the spatial average
            for cell in oneFire_cells:
                print('Working on cell {}...'.format(cell))
                data.loc[data['id'] == cell, 'season_q50s'] = spatial_average(cell, data, 'season_q50s', edge_list)
                data.loc[data['id'] == cell, 'season_q50s_rel'] = spatial_average(cell, data, 'season_q50s_rel', edge_list)
                
            
                
            # Update variable to contain entries with nans
            oneFire_cells = data.loc[np.isnan(data['season_q50s']), 'id']
            print(oneFire_cells)
                
                # # TEST 31/10/2023: what if we take the spatial average for the - DISCARDED
                # # other variables as well
                # for var in ['number_fires', 'burned_area', 'size_q50_int', 'size_q99_int', 'frequency_q99s']:
                #     # Create new column with same data
                #     data[var + '_sa'] = data[var].copy()
                    
                #     # Calculate the spatial average for cells with one fire
                #     data.loc[data['id'] == cell, var + '_sa'] = spatial_average(cell, data, var + '_sa', edge_list)
                
            
        # UPDATE 26/10/2023 Normalise by total flammable area the burned area and the 
        #                   number of fires
        data['burned_area_n'] = data['burned_area'] / data['flammable']
        data['number_fires_n'] = data['number_fires'] / data['flammable']
        
        # ---------------- Transforming variables ---------------------------------
        
        for i in range(len(list_vars_labels)):
            
            if list_vars_transform[i] == 'log':
                data['log_' + list_vars_labels[i]] = np.log(data[list_vars_labels[i]])
                
            elif list_vars_transform[i] == 'log1':
                data['log1_' + list_vars_labels[i]] = np.log(data[list_vars_labels[i]] + 1)
                
            elif list_vars_transform[i] == 'sqrt':
                data['sqrt_' + list_vars_labels[i]] = np.sqrt(data[list_vars_labels[i]])
                
            elif list_vars_transform[i] == 'log_rev':
                data['log_' + list_vars_labels[i]] = np.log(data[list_vars_labels[i]] + np.abs(np.min(data[list_vars_labels[i]])) + 1)
        

        # Saving file
        data.to_csv(out_fp.format(grid, year), index = False)
        


if __name__ == '__main__':
    main()
