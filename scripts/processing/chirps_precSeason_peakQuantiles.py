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
from readFormat_dataYears import readFormat_dataYears
from find_peakMonth import find_peak_perCell
from find_quantiles_inverseCDF import inverseCDF_perCell

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
        
        data_saved.to_csv(data_saved, index = False)
"""
MAIN
"""
def main():
    
    # ----------------- User inputs -------------------------------------------
    
    # Directory for input and output files
    indir_fp = '../../data/processed/climate_data/'
    out_dir = '../../data/processed/summary_tables/'
    
    # Files to monthly burned area per period
    data_fp = indir_fp + 'chirps_grid50km_{}.csv'
    
    # Output filepath
    output_fp = out_dir + 'precSeason_measures_chirps_period{}.csv'
    
    # Grid
    grid_id = '50km'
    # Minimum area to select cells
    cell_area_thresh = 75.0
    
    # Total number of periods we are working with
    nperiods = 4
    
    
    # --------------- Reading data --------------------------------------------
    
    print('Reading data...')
    data = readFormat_dataYears(data_fp, 1985, 2020, grid_id, ['mean_prec'], 
                             periods = True, nperiods = nperiods, 
                             thresh = cell_area_thresh).reset_index(drop = True)
    
    # --------------- Formatting data -----------------------------------------
    
    # Now we have total precipitation per month and year. We want average precipitation
    # per month and period
    print('Calculating average precpitation per month and period...')
    data = data[['id','%area_crop','period','month', 'mean_prec']].groupby(by = ['id','%area_crop','period','month']).mean().reset_index()

    
    # --------------- Getting the peak of the precipitaiton season ------------
    
    print('Finding the peak month per cell and period...')
    peak_season = find_peak_perCell(data, 'id', 'period', 'mean_prec', exclude_nan = False)
    
    
    # Subsetting columns of interest
    peak_season = peak_season[['id', 'period', 'peak_month']]
    

    
    # -------------------------------------------------------------------------
    
    
    
    # ------------- Getting the quantiles of the season -----------------------
    
    print('Finding the quantiles of the precipitation season...')
    quantiles = inverseCDF_perCell(data, 'id', 'period', 'mean_prec', exclude_nan = False)
    
    
    # -------------------------------------------------------------------------
    
    
    # ------------- Merging results -------------------------------------------
    
    print('Merging outputs...')
    output = peak_season.merge(quantiles, on = ['id','period'], how = 'outer')
    
    del data, peak_season, quantiles
    
    # -------------------------------------------------------------------------
    
    
    # # ------------------- Saving data by period -----------------------------
    
    print('Saving output...')    
    

    # Looping over the periods, saving the individual files
    for i, (period, data_subset) in list(zip(range(1, nperiods + 1), output.groupby('period'))):
        
        print('...for period {}'.format(period))
        save_fireSeasonCharacteristic(data_subset, output_fp.format(i))




if __name__ == '__main__':
    main()
