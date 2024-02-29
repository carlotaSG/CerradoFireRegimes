# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 20:04:26 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculates the fire season concentration index, and the month that 
cooresponds to the peak of the season according to this index.
"""

"""
Imports
"""
import pandas as pd
import numpy as np

from os.path import exists


import sys
sys.path.append('../utils/')
from readFormat_dataPeriods import readFormat_allPeriods
import burnedArea_concentrationIndex as baCI

"""
Functions
"""
def ba_concentrationIndex(ba_date, ba_col):
    # Function that calculate a burned area concentration index as an adaptation 
    # of the one used in Seasonality of Precipitation in the US and mentioned in Defining pyromes and ...
    
    # Burned Area Concentration INdex and peak of season according to this method
    ba_ci, peak_ci = baCI.burnedArea_concentrationIndex(ba_date, ba_col)
    
    # peak_ci is simply the month, get date instead
    # peak_ci = ba_date.loc[ba_date['month'] == peak_ci, 'date'].item()
    # print(ba_ci)
    # print(peak_ci)
    
    return ba_ci, peak_ci



def find_concentrationIndex_perCell(df, grouping_id, time_interval, var, exclude_nan = False):
    # This function receives a dataframe of the type:
    #   cell/grouping | time_interval | month | variable
    # gropus by cell/grouping and time_interval and calls ranking to find the 
    # legnth of the fire season - the number of months that add up to 80% of the burned area
    # Combinations of cell/grouping and year where variable is - for all months are
    # filled with NA and a message is produced
    
    # List to store the length (number of months) per cell and time_interval in the time series (ts)
    var_ci_ts = []
    
    # Looping over pair grouping - year
    for (grouping, time_int), group in df.groupby( by = [grouping_id, time_interval]):
        # It may be that there is no burned area in a certain period.
        # In this case, I will print an informative message and I will set the peak month as 0
        if group[var].sum() == 0:
            print('WARNING: no {} in cell {}, time interval {}'.format(var, grouping, time_int))
            ci, peak_ci = np.nan, np.nan
        else:
            ci, peak_ci = ba_concentrationIndex(group, var)
        
        # Appending peak month for this grouping and year to list of dataframes
        var_ci_ts.append(pd.DataFrame({
            grouping_id : grouping,
            time_interval : time_int,
            'ci' : [ci],
            'peak_ci': peak_ci
            }))
    
    # Converting list of dataframes to single dataframe
    var_ci_ts = pd.concat(var_ci_ts, ignore_index = True)
    
    # Exclude NaNs if asked to
    if exclude_nan == True:
        var_ci_ts = var_ci_ts[~ var_ci_ts['ci'].isna()].reset_index(drop = True)
        
    return var_ci_ts

def save_fireSeasonCharacteristic(data, data_fp):
    
    # First, check if output filepath already exists
    file_exists = exists(data_fp)
    
    # If it does not exist, save it directly
    if file_exists == False:
        data.to_csv(data_fp, index = False)
        
    # If it exists, read, merge data and save
    elif file_exists == True: 
        
        data_saved = pd.read_csv(data_fp)
        
        data_saved = data_saved.merge(data, on = ['id', 'period'], how = 'outer')
        
        
        data_saved.to_csv(data_fp, index = False)
"""
MAIN
"""
def main():
    
    # ----------------- User inputs -------------------------------------------
    
    # Number of years in period
    nyears = 9
    
    # Total number of periods we are working with
    nperiods = 4
    
    # Grid
    grid_id = '30km'
    
    # Directory for input and output files
    dir_fp = '../../data/processed/summary_tables/{}_year_periods/'.format(nyears)
    
    # Files to monthly burned area per period
    # data_fp = dir_fp + 'burnedArea_fireCount_intersection_monthFromPixel_period{}.csv'
    
    # # Output filepath
    # output_fp = dir_fp + 'fireSeason_measures_period{}.csv'
    
    data_fp = dir_fp + 'burnedArea_monthFromPixel_fireCount_assignedCell_exclude3RainyMonths_grid{}_period{}.csv'.format(grid_id, '{}')
    
    # Output filepath
    output_fp = dir_fp + 'fireSeason_measures_exclude3RainyMonths_grid{}_period{}.csv'
    
    
    # # Minimum area to select cells
    # cell_area_thresh = 75.0
    
    # # Extra cells to include and cells to exclude
    # include_cells = [567, 568]
    # exclude_cells = [1451]
    
    
    
    
    # --------------- Reading data --------------------------------------------
    
    print('Reading data...')
    data = readFormat_allPeriods(data_fp, nperiods, grid_id)#, thresh = cell_area_thresh, incl_cells = include_cells, excl_cells = exclude_cells)
    
    
    # --------------- Formatting data -----------------------------------------
    
    # To then identify the peak month
    # Converting the data from long to wide format
    # data = pd.melt(data, 
    #                id_vars = ['id','%area_crop', 'geometry','period','npolygons','area_T'], 
    #                var_name = 'month', value_name = 'area')
    
    # # Converting month column from string to integer
    # data['month'] = data['month'].apply(int)
    
    
    # --------------- Getting the peak of the fire season ---------------------
    
    print('Calculating the concentration index per cell and period...')
    concIndex = find_concentrationIndex_perCell(data, 'id', 'period', 'area', exclude_nan = False)
    
    
    # Subsetting columns of interest
    concIndex = concIndex[['id', 'period', 'ci', 'peak_ci']]
    
    # Freeing memory
    del data
    
    
    # ------------------- Saving data by period -------------------------------
    
    print('Saving output...')    
    
    # Cumulative variable indicating the period we are working on 
    i = 1
    # Looping over the periods, saving the individual files
    for period, data_subset in concIndex.groupby('period'):
        
        print('...for period {}'.format(period))
        
        # print(data_subset.head(20))
        save_fireSeasonCharacteristic(data_subset, output_fp.format(grid_id, i))
        
        i += 1





if __name__ == '__main__':
    main()
