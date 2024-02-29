# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 17:18:39 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculates the average precipitation per month, cell and period.

Input files:
    - ../../data/processed/climate_data/chirps_grid<grid_id>_<year>.csv
    
Output files:
    - ../../data/processed/summary_tables/<nperiods>_year_periods/chirps_meanMonthlyPrec_grid<>_period<>.csv
    

"""

"""
Imports
"""
import sys
sys.path.append('../utils/')

from readFormat_dataYears import readFormat_dataYears

"""
MAIN
"""
def main():
    
    # ----------------- User inputs -------------------------------------------
    
    # Total number of periods we are working with
    nperiods = 2
    nyears = 18
    
    # Grid id
    grid_id = '50km'
    
    start_year = 1985
    end_year= 2022
    
    # Directory for input and output files
    indir_fp = '../../data/processed/climate_data/'
    out_dir = '../../data/processed/summary_tables/{}_year_periods/'
    
    # Files to monthly burned area per period
    data_fp = indir_fp + 'chirps_grid{}_{}.csv'.format(grid_id, '{}')
    
    # Output filepath
    output_fp = out_dir + 'chirps_meanMonthlyPrec_grid{}_period{}.csv'.format(grid_id, '{}')
    

    # Minimum area to select cells
    cell_area_thresh = 0.0 # 75.0
    
    # # Extra cells to include and cells to exclude
    # include_cells = [567, 568]
    # exclude_cells = [1451]
    

    
    # --------------- Reading data --------------------------------------------
    
    print('Reading data...')
    data = readFormat_dataYears(data_fp, start_year, end_year, grid_id, ['mean_prec'], 
                             periods = True, nperiods = nperiods, 
                             thresh = cell_area_thresh).reset_index(drop = True)
    
    # --------------- Formatting data -----------------------------------------
    
    # Now we have total precipitation per month and year. We want average precipitation
    # per month and period
    print('Calculating average precpitation per month and period...')
    data = data[['id','%area_crop','period','month', 'mean_prec']].groupby(by = ['id','%area_crop','period','month']).mean().reset_index()
    
    
    # ---------------- Saving data --------------------------------------------
    
    print('Saving output...')
    
    for i, (period, data_subset) in list(zip(range(1, nperiods + 1), data.groupby('period'))):
        
        print('...for period {}'.format(period))
        data_subset.to_csv(output_fp.format(nyears, i), index = False)


if __name__ == '__main__':
    main()

