# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 11:20:12 2022

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that uses the information on the month with minimum average burned area
per grid cell in a time series to establish the sequence of fire years. 

A fire year is defined as in (REFERENCE), where the month with minimum burn area
is the last month in the fire year.

"""

"""
Imports
"""
import pandas as pd
import numpy as np

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
    
    # Period over which to work
    start_year = 1985
    end_year = 2020
    
    # UPDATE 15/2/2022
    # # File with burned area time series (month-year) per cell
    # timeseries_fp = '../data/processed/grid_{}/MBfogo_burnedArea_grid{}_cerrado_{}.csv'.format(grid_id, grid_id, period)
    
    
    # UPDATE 15/2/2022 
    # Changed to work with Precipitation too
    # File with minimum average burned area month per cell
    minBA_fp = '../data/processed/grid_{}/lookup_tables/lkup_month_minPrec_cell_{}-{}.csv'.format(grid_id, start_year, end_year)
    
    
    # ---------------- Reading data -------------------------------------------
    
    # Reading only time series columns - at this point we are onyl interested 
    # in having the full time series to which the information we have belongs
    # timeseries = pd.read_csv(timeseries_fp, usecols = ['year', 'month'])
    years = list(range(start_year, end_year + 1))
    months = list(range(1, 12+1))
    timeseries = pd.DataFrame({'year': [year for year in years for i in range(12)], 'month': months*len(years)})
    
    minBA = pd.read_csv(minBA_fp)
    
    # --------------------- Establishing fire years ---------------------------
    
    # List to accumulate dataframes with the information per cell
    list_df = []
    
    # For each cell in file (each row)
    for row in range(len(minBA)):
        
        # Create copy of timeseries dtaframe for this cell
        cell_timeseries = timeseries.copy()
        
        # Fire year -----------------------------------------------------------
        # First, create column where each year is relaballed as their position 
        # in the time series, starting at 0
        cell_timeseries['fire_year'] = cell_timeseries['year'] - cell_timeseries['year'].min()
        
        # Then, shift this information as many times as the month with minimum 
        # average area for this cell. This way, the fire years are adjusted to 
        # finish in the month with minimum average burned area.
        cell_timeseries['fire_year'] = cell_timeseries['fire_year'].shift(periods = minBA.loc[row, 'month_minBA'])
        
        # Fire month ----------------------------------------------------------
        # Create column for fire year month, as the month column shifted 
        # as many times as the month with minimum average area for this cell.
        cell_timeseries['fire_month'] = cell_timeseries['month'].shift(periods = minBA.loc[row, 'month_minBA'])
        
        
        # Cell ID -------------------------------------------------------------
        # Create column for this cell
        cell_timeseries['id'] = minBA.loc[row, 'id']
        
        # Append dataframe to list
        list_df.append(cell_timeseries)
 
    # Convert list of dtaframes to dataframes
    timeseries = pd.concat(list_df, ignore_index = True)
    
    # Save to file
    timeseries.to_csv('../data/processed/grid_{}/lookup_tables/lkup_prec-year_timeseries_cell_{}-{}.csv'.format(grid_id, start_year, end_year), index = False)        
        
                            
    
if __name__ == '__main__':
    main()

