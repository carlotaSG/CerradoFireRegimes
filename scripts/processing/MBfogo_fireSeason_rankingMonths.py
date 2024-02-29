# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 20:04:26 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculates the length of the fire season following (Archibald et al. 2013)
ranking months by burned area and gettin ght enumber of month containing x% of 
the period's burned area.
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

"""
Functions
"""
def ranking(ba_date, ba_col, ranking_thresh):
    # Function that calculates the length of the fire season by ranking months in the 
    # year and counting how many months are needed to add up to 80% of total
    # burned area in the year.

    # First, sort dataframe by burned area values from highest to lowest
    ba_date_sorted = ba_date.sort_values(by = ba_col, ascending = False).reset_index(drop = True)
    # print(ba_date_sorted)
    # Then, calculate inverse cumulative percentage
    cumPerc = 100 * (ba_date_sorted[ba_col].cumsum() / ba_date_sorted[ba_col].sum())
    # print(cumPerc)
    # Get season length as the length of the subset That first adds up to 80%
    # We choose as the cutoff month, the month closest to 80% cumulatively.
    the_index = cumPerc.sub(ranking_thresh).abs().idxmin()
    season_length = the_index + 1
    # print(the_index)
    # If the difference between the percentage burned area for this month and 80%
    # is larger than 5%, then we choose the first month with percentage smaller than 80%
    if cumPerc[the_index] - ranking_thresh > 2:
        the_index = the_index - 1
        
        if the_index < 0:
            season_length = 1
            
        elif the_index >= 0:
            season_length = the_index + 1
    # If there is no such month, then we choose the first month as the only month in the
    # burn season, hence length = 1
    # print(season_length)
    
    return season_length




def find_length_perCell(df, grouping_id, time_interval, var, ranking_thresh, exclude_nan = False):
    # This function receives a dataframe of the type:
    #   cell/grouping | time_interval | month | variable
    # gropus by cell/grouping and time_interval and calls ranking to find the 
    # legnth of the fire season - the number of months that add up to 80% of the burned area
    # Combinations of cell/grouping and year where variable is - for all months are
    # filled with NA and a message is produced
    
    # List to store the length (number of months) per cell and time_interval in the time series (ts)
    var_length_ts = []
    
    # Looping over pair grouping - year
    for (grouping, time_int), group in df.groupby( by = [grouping_id, time_interval]):
        # It may be that there is no burned area in a certain period.
        # In this case, I will print an informative message and I will set the peak month as 0
        if group[var].sum() == 0:
            print('WARNING: no {} in cell {}, time interval {}'.format(var, grouping, time_int))
            length_season = np.nan
        else:
            length_season = ranking(group, var, ranking_thresh)
        
        # Appending peak month for this grouping and year to list of dataframes
        var_length_ts.append(pd.DataFrame({
            grouping_id : grouping,
            time_interval : time_int,
            'length_ranking_{}'.format(ranking_thresh) : [length_season]
            }))
    
    # Converting list of dataframes to single dataframe
    var_length_ts = pd.concat(var_length_ts, ignore_index = True)
    
    # Exclude NaNs is asked to
    if exclude_nan == True:
        var_length_ts = var_length_ts[~ var_length_ts['peak_month'].isna()].reset_index(drop = True)
        
    return var_length_ts

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
    
    # Grid
    grid_id = '30km'
    
    # Years in period
    nyears = 9
    
    # Total number of periods we are working with
    nperiods = 4
    
    # Directory for input and output files
    dir_fp = '../../data/processed/summary_tables/{}_year_periods/'.format(nyears)
    
    # Files to monthly burned area per period
    # data_fp = dir_fp + 'burnedArea_fireCount_intersection_monthFromPixel_period{}.csv'
    
    # # Output filepath
    # output_fp = dir_fp + 'fireSeason_measures_period{}.csv'
    
    data_fp = dir_fp + 'burnedArea_monthFromPixel_fireCount_assignedCell_exclude3RainyMonths_grid{}_period{}.csv'.format(grid_id, '{}')
    
    # Output filepath
    output_fp = dir_fp + 'fireSeason_measures_exclude3RainyMonths_grid{}_period{}.csv'
    
    
    # Minimum area to select cells
    # cell_area_thresh = 75.0
    
    
    
    # The threshold % to consider the length of the fire seaons in the ranking 
    # (that is, the number of months that contain the ranking_threshold)
    ranking_threshold = 80
    
    
    # --------------- Reading data --------------------------------------------
    
    print('Reading data...')
    data = readFormat_allPeriods(data_fp, nperiods, grid_id)#, thresh = cell_area_thresh)
    
    
    # --------------- Formatting data -----------------------------------------
    
    # To then identify the peak month
    # Converting the data from long to wide format
    # data = pd.melt(data, 
    #                id_vars = ['id','%area_crop', 'geometry','period','npolygons','area_T'], 
    #                var_name = 'month', value_name = 'area')
    
    
    # --------------- Getting the peak of the fire season ---------------------
    
    print('Finding the length of the season per cell and period...')
    length_season = find_length_perCell(data, 'id', 'period', 'area', ranking_threshold, exclude_nan = False)
    
    
    # Subsetting columns of interest
    length_season = length_season[['id', 'period', 'length_ranking_{}'.format(ranking_threshold)]]
    
    # Freeing memory
    del data
    
    
    # ------------------- Saving data by period -------------------------------
    
    print('Saving output...')    
    
    # Cumulative variable indicating the period we are working on 
    i = 1
    # Looping over the periods, saving the individual files
    for period, data_subset in length_season.groupby('period'):
        
        print('...for period {}'.format(period))
        
        # print(data_subset.head(20))
        save_fireSeasonCharacteristic(data_subset, output_fp.format(grid_id, i))
        
        i += 1





if __name__ == '__main__':
    main()
