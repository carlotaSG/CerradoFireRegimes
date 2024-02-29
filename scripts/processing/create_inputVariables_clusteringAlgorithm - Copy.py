# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 10:19:31 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

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
    - concentration index (of the burned area season)
    - burned area (of the polygons intersection)
    
Input files:
    
    
Output file:
    ../../data/processed/clustering_algorithm/input_vars_grid{}_period{}.csv
    

"""

"""
Imports
"""
import pandas as pd
from functools import reduce


"""
Functions
"""
def replace(my_list, my_dict):
    return [x if x not in my_dict else my_dict[x] for x in my_list]


"""
MAIN
"""
def main():
    
    # -------------------------- User inputs ----------------------------------
    
    # Grid id
    grid = '50km'
    
    # Periods
    periods = list(range(1, 4+1))
    
    # Files
    #--------
    
    # Input directory
    in_dir = '../../data/processed/summary_tables/'
    
    # Filename end
    file_termination = '_grid{}_period{}.csv'
    
    # File identifiers
    fire_size_fp = 'fireSize_assigned'                                          # Right subset
    number_fires_fp = 'burnedArea_fireCount_assigned'                           # All cells
    frequency_fp = 'fireFrequency_fromPixels_sECDF_quantiles'                   # 'fireFrequency_fromPixels_quantiles'        
    fire_return_interval_fp = 'fireFrequency_returnInterval'                 
    season_fp = 'fireSeason_measures_exclude3RainyMonths'
    burned_area_fp = 'burnedArea_fireCount_intersection_monthFromPolygon'
    
    # Cast file identifiers into list
    file_ids = [fire_size_fp, number_fires_fp, frequency_fp, fire_return_interval_fp, season_fp, burned_area_fp]
    
    
    # Output filepath
    out_fp = '../../data/processed/clustering_algorithm/input_vars_sECDF_grid{}_period{}.csv'
    
    
    # Variables
    #-----------
    fire_season_var = ['ci', 'median_month', 'median_month_sECDF'] # ['ci', 'median_month']
    fire_size = ['q_0.5', 'q_0.99', 'pwl_alpha']
    number_fires_var = ['npolygons']
    burned_area_var = ['area_T']
    frequency_var = ['f99']
    fri_var = ['fri']
    
    # Cast variables into list
    list_vars = fire_season_var + fire_size + number_fires_var + burned_area_var + frequency_var + fri_var
    
    
    # Columns names
    #---------------
    
    # Adjusting column names to desired ones
    new_col_names = {'ci': 'season_ci',
                     'median_month': 'season_q50',
                     'median_month_sECDF': 'season_q50s',
                     'q_0.5': 'size_q50',
                     'q_0.99': 'size_q99',
                     'pwl_alpha': 'size_alpha',
                     'npolygons': 'number_fires',
                     'area_T': 'burned_area',
                     'f99': 'frequency_q99'}
    
    
    # Assigning mean neighbour value to cell 240
    c240_neighbours = [186, 241, 294]
    
    # ------------------ Reading and formatting input data --------------------
    
    # Working on one period at a time
    for period in periods:
        
        print('Working on period {}'.format(period))
        
        data = []
        
        # Reading files consecutively
        for file_fp in file_ids:
            data_file = pd.read_csv(in_dir + 
                                    file_fp + 
                                    file_termination.format(grid, period))
            
            # Manual cleaning to avoid same-column names. 
            if file_fp == number_fires_fp:
                data_file = data_file.drop(columns = 'area_T')
                
            elif file_fp == burned_area_fp:
                data_file = data_file.drop(columns = 'npolygons')
            
            elif file_fp == fire_return_interval_fp:
                data_file = data_file.drop(columns = ['area_T', 'flammable'])
            
            # Selecting relevant columns only
            data_file = data_file[data_file.columns.intersection(['id', 'period'] + list_vars)]
            
            # Adjusting names of columns
            data_file.columns = replace(data_file.columns, new_col_names)
            
            # Appending to list of dataframes
            data.append(data_file)
            
            # Freeing memory
            del data_file
        
        # Casting list of dataframes to unique dataframe
        data = reduce(lambda x, y: x.merge(y, on = ['id','period']), data)
        
        # Creating column with the frequency from fire return interval
        data['fri_freq'] = data['fri'].apply(lambda x: 9/x)
        
        # Function reduce(function, iterable) applies function of two arguments
        # cumulatively to the items of the iterable from left to right
        
        if period == 1:
            # Before saving, cell 240 in period 1 does not have any measures associated, 
            # so we are going to "fill it" with the average of the cells around it
            
            # First, convert column id to index
            data = data.set_index('id')
            
            # Then, add row with index 240
            data.loc[240] = ['(1985, 1993)'] + \
                data.loc[data.index.isin(c240_neighbours), data.columns != 'period'].apply('mean', axis = 0).tolist()
            
            # Reset index
            data = data.reset_index()
            
            # Sort according to column id
            data = data.sort_values(by = 'id').reset_index(drop = True)
            
            data.loc[data['id'] == 240, 'number_fires'] = int(data.loc[data['id'] == 240, 'number_fires'])
            data.loc[data['id'] == 240, 'season_q50'] = int(data.loc[data['id'] == 240, 'season_q50'])
        

        # Saving file
        data.to_csv(out_fp.format(grid, period), index = False)
        


if __name__ == '__main__':
    main()
