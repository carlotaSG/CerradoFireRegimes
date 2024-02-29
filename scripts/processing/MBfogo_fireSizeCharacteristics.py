# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 16:05:11 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

**IMPORTANT: This script must be run in the scipy-env environment.

Script that calculates the fire size characteristics per cell and period. The 
fire characteristics calculated are:
    - 5% quantile
    - 50% quantile (median)
    - 95% quantile
    - the exponent (alpha) and sigma (standard error) of a fitted power law
    
    
Input files:
    - ../../data/processed/MBfogo_c10_cerrado_polygons_b150_CM/MBfogo_c10_cerrado_polygons_{}.csv

Output files:
    - ../../data/processed/summary_tables/fireSize_assigned_grid{}_period{}.csv'   

"""

"""
Imports
"""
import pandas as pd

import powerlaw

import sys
sys.path.append('../utils/')
from readFormat_dataYears import readFormat_dataYears


"""
Functions
"""
def getQuantiles(df, list_q):
    """
    Function that calculates the quantiles of fire size (area) indicated in the 
    list of quantiles, list_q, for each cell in dataframe. Returns a dataframe
    of the form:
        id | q_<>

    Parameters
    ----------
    df : DataFrame
        Dataframe containing a list of fire  polygons, each one with the grid cell
        they belong to (id) and their aream in km2 (area).
    list_q : List of floats
        The different quantiles we want to describe the fire size distribution in each cell.

    Returns
    -------
    df_cell_quantiles : DataFrame
        Dataframe of the form: id | q_<>

    """
    
    # Grouping by cell and getting the different quantiles
    df_cell_quantiles = df[['id', 'area']].groupby(by = 'id').quantile(q = list_q).reset_index()
    
    # Casting from wide to long
    df_cell_quantiles = df_cell_quantiles.pivot(index = 'id', columns = 'level_1', values = 'area').reset_index()
    
    df_cell_quantiles.columns = ['id'] + ['q_'+ str(q) for q in list_q]
    
    return df_cell_quantiles


def fitPowerlaw(df, x_min):
    """
    Function that fits a powerlaw to a fire size distribution. Returns teh input
    dataframe with two additional columns indicating the the exponent of the
    distribution and its standard error.

    Parameters
    ----------
    df : DataFrame
        Dataframe where each row is a different polygon and contains the columns
        'id' indicating the cell the polygon belongs to, and its 'area' in km2.
    x_min : float
        The minimu possible size of the fire polygon size distribution.

    Returns
    -------
    df : DataFrame
        Dataframe where each row is a cell, with the fitted power-law exponent
        alpha and its standard error. Of the form: id | pwl_alpha | pwl_sigma.

    """
    
    # Group by cell and create powerlaw.Fit object for each group of areas, 
    # passing the xmin value provided by user
    df = df[['id', 'area']].groupby(by = 'id').agg({'area': lambda x: powerlaw.Fit(x, xmin = x_min)})
    
    # Changing column name for clarity
    df = df.rename(columns = {'area': 'fit'})
    
    # Calculating the exponent (alpha) and its standard error (sigma)
    df['pwl_alpha'] = df['fit'].apply(lambda x: x.power_law.alpha)
    df['pwl_sigma'] = df['fit'].apply(lambda x: x.power_law.sigma)
    df = df.drop(columns = 'fit').reset_index()
    
    return df

"""
MAIN
"""
def main():
    
    # ---------------- User inputs --------------------------------------------
    
    # Grid to work on 
    grid_id = '30km'
    grid_col_name = 'cellID_{}'.format(grid_id)
    
    # Number of periods
    n_periods = 4
    
    if n_periods == 2:
        # Year periods to work on 
        periods = [(1985, 2002), (2003, 2020)]
        
        # Number of years per period
        n_years = 18
        
    elif n_periods == 4:
        # Year periods to work on 
        periods = [(1985,1993), (1994, 2002), (2003, 2011), (2012, 2020)]
        
        # Number of years per period
        n_years = 9
    
    period_labels = [str(p) for p in periods]
    
    # Start and end years
    start_year = 1985
    end_year = 2020
    
    # UPDATE (24/10/2023): Allowing to choose the MapBiomas FOGO collection
    collection = 2
    
    # Filepath to files
    data_year_fp = '../../data/processed/MBfogo_c{}0_cerrado_polygons_b150_CM/MBfogo_c{}0_cerrado_polygons_{}.csv'.format(collection, collection, '{}')
    
    # Output file
    output_fp = '../../data/processed/summary_tables/{}_year_periods_newGrids/fireSize_assigned_grid{}_period{}.csv'


    # UPDATE (24/10/2023): We calculate the data for all grid cells
    # Work only with cell that have more than a certain percentage of their area inside the Cerrado
    thresh = 0.0 # 75.0
    
    # # Extra cells to include and cells to exclude
    # include_cells = [567, 568]
    # exclude_cells = [1451]
    
    # List of quantiles we want to work on
    list_quantiles = [0.05, 0.5, 0.95, 0.99, 1.]
    
    # The minimum polygon area that we accept
    min_polSize = 0.03
    
    
    
    #--------------------------------------------------------------------------
    
    # ------------- Reading data ----------------------------------------------
    
    print('Reading data for all years...\n')
    
    # Read and format the data for all years consecutively
    data_year = readFormat_dataYears(data_year_fp, start_year, end_year, grid_id, ['area'], 
                                     periods = True, nperiods = n_periods, grid_colname = grid_col_name, 
                                     fire_polygons = True, min_size = min_polSize, 
                                     thresh = thresh)
    
        
    # -------------------------------------------------------------------------
    
    # -------------- Fire characteristics -------------------------------------
    
    i = 0
    # Calculating the fire size characteristics for each cell and period
    for period in data_year['period'].unique():
        
        print('Working on period {}...'.format(period))
        
        # Subsetting the data for this period
        data_period = data_year[data_year['period'] == period]
        
        
        # Calculating area quantiles for each cell
        data_quantiles = getQuantiles(data_period, list_quantiles)
        
        # Fitting a powerlaw to each cell's areas and calculating the exponent
        # the corresponding standard error
        data_pwl = fitPowerlaw(data_period, min_polSize)
        
        # Merging the two bits of data
        data_period = data_quantiles.merge(data_pwl, on = 'id', how = 'outer')
        
        # Freeing memory
        del data_quantiles, data_pwl
        
        # Adding column indicating the period
        data_period['period'] = period_labels[i]
        i += 1
        
        # Reordering columns
        data_period = data_period[['id', 'period'] + data_period.columns[1:-1].tolist()]
        
        # Saving to file
        data_period.to_csv(output_fp.format(n_years, grid_id, i), index = False)
        
        print('')
    
    
    
    
if __name__ == '__main__':
    main()
