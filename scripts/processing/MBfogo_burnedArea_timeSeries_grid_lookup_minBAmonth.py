# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 10:34:55 2022

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Short script that gets output from MBfogo_burnedArea_timeSeries_grid.py and 
calculates, for each cell, the month of minimum average activity in the time series.
Saves information to an output look-up table.

"""

"""
Imports
"""
import pandas as pd

"""
Functions
"""

"""
MAIN
"""
def main():
    
    # -------------------- User inputs ----------------------------------------
    
    # Grid size
    grid_size = 100
    # Units of grid_size
    grid_units = 'km'
    
    grid_id = str(grid_size) + grid_units
    
    # # Period over which to work
    # period = '1990-2020'
    
    # # File with burned area time series (month-year) per cell
    # data_fp = '../data/processed/grid_{}/MBfogo_burnedArea_grid{}_cerrado_{}.csv'.format(grid_id, grid_id, period)
    
    # UPDATE 15/02/2022
    start_year = 1985
    end_year = 2020
    
    data_fp = '../data/processed/grid_{}/MBfogo_burnedArea_grid{}_cerrado_{}.csv'
    
    # Update 15/2/2022
    # Doing it also for precipitation
    data_fp = '../data/processed/grid_{}/summary_tables/summary_era5land_{}_{}.csv'
    
    var_col = 'mean_prec'
    
    
    # ----------------------- Calculations ------------------------------------
    
    # # Read data
    # data = pd.read_csv(data_fp)
    
    # Read data
    data = []
    for year in range(start_year, end_year + 1):
        data_year = pd.read_csv(data_fp.format(grid_id, grid_id, year))
        
        if 'year' not in data_year.columns:
            data_year['year'] = year
        
        data.append(data_year)
    data = pd.concat(data, ignore_index = True)
    
    
    
    # Reshaping from long to wide
    data = data.pivot(index = ['year', 'month'], columns = 'id', values = var_col).reset_index()
    data = data.dropna(axis = 1, how ='all')
    # print(data)
    # Calculate time-series average burned area per month
    # First, group by month and calculate average burned area (dropping year column)
    data = data.iloc[:, 1:].groupby(by = 'month').mean().reset_index()
    
    # Identifying the month of minimum average burned area per cell
    data = data.iloc[:, 1:].apply(lambda x: data.loc[x.idxmin(), 'month']).reset_index()
    data.columns = ['id', 'month_minBA']
    
    # Saving information to lookup table
    data.to_csv('../data/processed/grid_{}/lookup_tables/lkup_month_minPrec_cell_{}-{}.csv'.format(grid_id, start_year, end_year), index = False)
    
    
    
    
    
if __name__ == '__main__':
    main()
