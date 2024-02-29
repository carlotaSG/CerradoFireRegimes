# -*- coding: utf-8 -*-
"""
Created on Sun Feb 27 18:51:00 2022

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that relabels months in a dataframe column from a month of reference to 
which a new label is given. This script is used to relabel months so that a linear
regression can be adjusted to this circular data making sure January goes after
December and vice versa.

It does so working on grid cells consecutively.
"""

"""
Imports
"""
import pandas as pd
import geopandas as gpd

"""
Functions
"""
def adjustMonth(month):
    """
    Once the initial set of months are transformed by adding the transformation 
    value, some months may be negative and others may be larger than 12. Hence, 
    these months need to be corrected back to being between 1 and 12. This is 
    done considering that January goes after December, that is, that a fire 
    year is a cycle.

    Parameters
    ----------
    month : list (Pandas series)
        List of transformed months that need to be corrected.
    transformation : integer
        Tnumber by which the months have been transformed.

    Returns
    -------
    month : list (Pandas series)
        List of corrected transformed months. Values between 1 and 12.

    """
    
    if month > 12:
        month = month - 12
    elif month < 1:
        month = month + 12
        
    return month

"""
MAIN
"""
def main():
    
    # User inputs -------------------------------------------------------------
    
    # Path to file with list of months to be used as reference for each month
    ref_fp = '../data/processed/grid_100km/lookup_tables/lkup_month_minmax_ba_km2_cell_1985-2020_MBfogo.csv'
    
    # Variable containing the month reference
    ref_var = 'month_max_ba_km2'
    
    # Path to file with the month time series
    month_ts_fp = '../data/processed/grid_100km/fire-season_definition/fire-season_data_grid100km_cell_1985-2020.csv'
    
    # Variable containing the month time series
    month_ts_var = 'peak_season'
    
    # Grid data filepath - to be used to select cells with more than 70% of their area within Cerrado
    grid_fp = '../data/grids/cerrado_grid_100km_cropped.shp'
    
    # Output filepath
    output_fp = '../data/processed/grid_100km/fire-season_definition/fire-season_data_grid100km_cell_1985-2020_peakSeason_ts_relabelled.csv'
    
    # Reading files -----------------------------------------------------------
    
    # Grid
    grid = gpd.read_file(grid_fp)
    list_cells = grid.loc[grid['%area_crop'] >= 70., 'id'].tolist()
    
    # Reference month per cell
    ref = pd.read_csv(ref_fp)
    # Filtering cells
    ref = ref[ref['id'].isin(list_cells)]
    # print(ref)
    
    # Month time series
    month_ts = pd.read_csv(month_ts_fp)
    month_ts = month_ts[month_ts['id'].isin(list_cells)]
    
    # Selecting and formatting the time series
    month_ts = month_ts[['id', 'fyear', month_ts_var]]
    
    month_ts[month_ts_var] = pd.to_datetime(month_ts['peak_season'], format =  '%d/%m/%Y')
    
    month_ts['peak_season_month'] = pd.DatetimeIndex(month_ts[month_ts_var]).month
    month_ts['year'] = pd.DatetimeIndex(month_ts[month_ts_var]).year
    
    # print(month_ts)
    
    # List of datafames with relabelled time series for each cell
    month_ts_relabel = []
    
    # For each cell, calculate how much to offset so that peak season month is 6
    # Add this number of the time-series of months
    for cell in month_ts['id'].unique():
        
        print('Working on cell {}'.format(cell))
        
        # Getting the average peak of the BA season for this cell
        peak_month_cell = ref.loc[ref['id'] == cell, ref_var].item()
        
        # Calculate how much do we have to add to this value so that this month becomes 6
        offset = 6 - peak_month_cell
        
        # Subset time series for this cell
        month_ts_cell = month_ts[month_ts['id'] == cell].reset_index(drop = True)
        
        # Fill gaps with average if there are NaNs
        if month_ts_cell['year'].isnull().any():
            
            # List of indices where there is a NaN
            nan_i = month_ts_cell[month_ts_cell['year'].isnull()].index
            
            for i in nan_i:
                if i != 35:
                    month_ts_cell.loc[i, 'year'] = int((month_ts_cell.loc[i+1, 'year'] + month_ts_cell.loc[i-1, 'year'])/2)
                    month_ts_cell.loc[i, 'peak_season_month'] = int((month_ts_cell.loc[i+1, 'peak_season_month'] + month_ts_cell.loc[i-1, 'peak_season_month'])/2)
                if i == 35:
                    month_ts_cell.loc[i, 'year'] = int(month_ts_cell.loc[i-1, 'year'])
                    month_ts_cell.loc[i, 'peak_season_month'] = int(month_ts_cell.loc[i-11, 'peak_season_month'])
        
        # Adding offset to months
        month_ts_cell['peak_season_month_relabel'] = month_ts_cell['peak_season_month'] + offset
        
        # print(month_ts_cell)
        
        # Correcting months that now are negative or larger then 12
        month_ts_cell['peak_season_month_relabel'] = month_ts_cell['peak_season_month_relabel'].apply(adjustMonth)
        
        # print(month_ts_cell)
        
        # Add this dataframe to list
        month_ts_relabel.append(month_ts_cell)
        
        print("")
        # break
    month_ts_relabel = pd.concat(month_ts_relabel, ignore_index = True)
    
    month_ts_relabel.to_csv(output_fp, index = False)

if __name__ == '__main__':
    main()