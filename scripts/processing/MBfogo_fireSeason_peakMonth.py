# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 20:04:26 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that identifies the peak of the fire season (the month with largest 
burned area) for each cell and period. Saves the output into one file per period.
"""

"""
Imports
"""
import pandas as pd

from os.path import exists


import sys
sys.path.append('../utils/')
from readFormat_dataPeriods import readFormat_allPeriods
from find_peakMonth import find_peak_perCell

"""
Functions
"""
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
    
    # Total number of periods we are working with
    nperiods = 4
    
    # Number of years in period
    nyears = 9
    
    # Directory for input and output files
    dir_fp = '../../data/processed/summary_tables/{}_year_periods/'.format(nyears)
    
    # Files to monthly burned area per period
    data_fp = dir_fp + 'burnedArea_monthFromPixel_fireCount_assignedCell_exclude3RainyMonths_grid{}_period{}.csv'.format(grid_id, '{}')
    # data_fp = dir_fp + 'burnedArea_fireCount_intersection_monthFromPixel_period{}.csv'
    
    # Output filepath
    output_fp = dir_fp + 'fireSeason_measures_exclude3RainyMonths_grid{}_period{}.csv'
    
    
    # # Minimum area to select cells
    # cell_area_thresh = 75.0
    
    # # Extra cells to include and cells to exclude
    # include_cells = [567, 568]
    # exclude_cells = [1451]
    
    # Total number of periods we are working with
    nperiods = 4
    
    
    # --------------- Reading data --------------------------------------------
    
    print('Reading data...')
    data = readFormat_allPeriods(data_fp, nperiods, grid_id)#, thresh = cell_area_thresh, incl_cells = include_cells, excl_cells = exclude_cells)
    
    
    # --------------- Formatting data -----------------------------------------
    
    # To then identify the peak month
    # Converting the data from long to wide format
    # data = pd.melt(data, 
    #                id_vars = ['id','%area_crop', 'geometry','period','npolygons','area_T'], 
    #                var_name = 'month', value_name = 'area')
    
    
    # --------------- Getting the peak of the fire season ---------------------
    
    print('Finding the peak month per cell and period...')
    peak_season = find_peak_perCell(data, 'id', 'period', 'area', exclude_nan = False)
    
    
    # Subsetting columns of interest
    peak_season = peak_season[['id', 'period', 'peak_month']]
    
    # Freeing memory
    del data
    
    
    # ------------------- Saving data by period -------------------------------
    
    print('Saving output...')    
    
    # Cumulative variable indicating the period we are working on 
    i = 1
    # Looping over the periods, saving the individual files
    for period, data_subset in peak_season.groupby('period'):
        
        print('...for period {}'.format(period))
        
        # print(data_subset.head(20))
        save_fireSeasonCharacteristic(data_subset, output_fp.format(grid_id, i))
        
        i += 1



if __name__ == '__main__':
    main()
