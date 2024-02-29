# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 13:13:32 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculates the fire season characteristics 5th quantile, median
and 95th quantile using the smoothed empirical cumulative distribution 
function (ECDF).


This program is very similar to the MBfogo_fireSeason_cumulativeProbability.py
but in this case it returns a continuous variable instead of discrete.

Input files:
    - ../../data/processed/summary_tables/{}_year_periods/burnedArea_monthFromPixel_fireCount_assignedCell_exclude3RainyMonths_period<>.csv
    
    
Output files:
    - ../../data/processed/summary_tables/{}_year_periods/fireSeason_measures_exclude3RainyMonths_grid<>_period<>.csv

"""

"""
Imports
"""
import pandas as pd
import numpy as np
import statsmodels.distributions.empirical_distribution as edf
from scipy.interpolate import interp1d

from os.path import exists

import sys
sys.path.append('../utils/')
from readFormat_dataPeriods import readFormat_allPeriods

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

def relabelMonths(df, var):
    """
    Function that relabels the months so that the month with largest burned area
    is in the middle of the fire year (it gets assigned label 6). The rest of 
    months get relabelled accordingly. Returns a dataframe with the initial and
    the new month labels, as well as the offset (number of month to be added to
    relabel the months.)

    Parameters
    ----------
    df : DataFrame
        Dataframe containing a column var and a month. The column month is the 
        one that is going to be relabelled (by creating a new column), and the
        var column is the one to use as reference to determine which month is to 
        be considered as month 6 (middle of the fire year).
    var : string
        The column to take as reference. The month corresponding to the max var
        value will be assigned number 6.

    Returns
    -------
    df : DataFrame
        Same dataframe as input with an additional column containing the relabelled
        months.
    offset : integer
        The number of months to add or subtract from the original months to convert
        original months to relabelled and vice versa.
    """
    
    # Calculate month with largest burned area
    peak_month = df.loc[df[var].idxmax(), 'month']
    
    # Calculate how much do we have to add to this value so that this month becomes 6
    offset = 6 - peak_month
    
    # Adding offset to months (relabelling)
    df['month_relabel'] = df['month'] + offset
    
    # Correcting months that now are negative or larger then 12
    df['month_relabel'] = df['month_relabel'].apply(adjustMonth)
    
    return df, offset



def calculate_sECDF(data):
    """
    Function that calculate the smoothed empirical cumulative distribution (sECDF)
    function from a list of values whose frequency describes the density function.
    
    It calcualtes the ECDF using the library statsmodels.distributions.empirical_distribution

    Parameters
    ----------
    data : list
        List of values from which to calcualte the sECDF.

    Returns
    -------
    q05 : float
        The 5th quantile of the sECDF.
    q50 : float
        The 50th quantile of the sECDF.
    q95 : float
        The 95th quantile of the sECDF.
    length : float
        The difference between q95 and 95, which would be the length of the 
        fire season in months.

    """
    
    # First, fit the ecdf
    data_ecdf = edf.ECDF(data)
    
    slope_changes = sorted(set(data))
    slope_changes = list(range(1, 12 + 1))
    
    data_ecdf_values_at_slope_changes = [data_ecdf(item) for item in slope_changes]
    
    inverted_ecdf = interp1d(data_ecdf_values_at_slope_changes,
                             slope_changes)
    
    # It may be that the 5th quantile is below the interpolation range's minimum
    # value. In that case, we use this value instead of the 5th quantile
    if 0.05 < data_ecdf_values_at_slope_changes[0]:
        q05val = data_ecdf_values_at_slope_changes[0]
    else:
        q05val = 0.05
    
    q05 = inverted_ecdf(q05val)
    q50 = inverted_ecdf(0.5)
    q95 = inverted_ecdf(0.95)
    
    length = q95 - q05
    
    return q05, q50, q95, length

def calculate_monthDistribution(data, col_var):
    """
    Function that gets data containing the burned area per month and cell.
    It converts the area to number of pixels and then generates a list of relabelled
    months, where each month appears as many times as the number of pixels that 
    burned in it. Then, it calls the function calculate_sECDF to calculate
    the quantiles of the smoothed Empirical Cumulative Distribution Function and
    the length of the fire season. It returns these values in the original-months
    units.

    Parameters
    ----------
    data : DataFrame
        Dataframe containing the burned area per month and cell.
    col_var : string
        The name of the column to use as a reference to relabel the months
        so that the distribution of col_var is cenetered on the month with maximum var.

    Returns
    -------
    m05 : float
        The 5th quantile of the sECDF.
    m50 : float
        The 50th quantile of the sECDF.
    m95 : float
        The 95th quantile of the sECDF.
    length : float
        The difference between q95 and 95, which would be the length of the 
        fire season in months.

    """
    
    # Make sure data is ordered by month
    data = data.sort_values(by = 'month').reset_index(drop = True)
    
    # Relabel the months so that peak month is number 6 (distribution is 
    # centered around the peak and in the middle of the fire year)
    data, offset = relabelMonths(data, col_var)
    
    # Convert are to number of pixels
    data['npixels'] = data['area']*1e06/900
    data['npixels'] = data['npixels'].apply(lambda x: int(x))
    
    # Now we create a vector of relabelled months containing each month repeated
    # as many times as the pixel
    relabel_month_data = []
    for rm in data['month_relabel'].unique():
        
        relabel_month_data = relabel_month_data + \
            [rm] * data.loc[data['month_relabel'] == rm, 'npixels'].item()
    
    
    # Calculate the quantiles obtained when fitting a smoothed ECDF,
    # and the length of the fire season
    m05, m50, m95, length = calculate_sECDF(relabel_month_data)
    
    
    # 4. Convert these quantiles back to the original months
    m05 = adjustMonth(m05 - offset)
    m50 = adjustMonth(m50 - offset)
    m95 = adjustMonth(m95 - offset)
    
    
    return m05, m50, m95, length

def inverse_smoothedECDF_perCell(df, grouping_id, time_interval, var, exclude_nan = False):
    """
    This function receives a dataframe in the format:
        cell/grouping | time_interval | month | quantiles
    groups by cell/grouping and time_interval and calls functions to calculate
    the quantiles of the fire season using the smoothed Empirical Cumulative 
    Distribution Function centered about the month with maximum area (var).

    Parameters
    ----------
    df : DataFrmae
        Dataframe containing the burned area per month, period and cell.
    grouping_id : string
        The name of the column by which to group the data (cell id).
    time_interval : string
        The name of the column containing the period identifier by which to group
        the data.
    var : string
        The column to use as reference to center the fire season on (e.g. the 
        burned area).
    exclude_nan : boolean, optional
        Whether to exclude the NaN values. The default is False.

    Returns
    -------
    var_quant_ts : DataFrame
        Dataframe containing the quantiles and length of the fire monthly 
        distribution per grouping_id and time_interval.

    """

    # List to store the length (number of months) per cell and time_interval in the time series (ts)
    var_quant_ts = []
    
    # Looping over pair grouping - year
    for (grouping, time_int), group in df.groupby( by = [grouping_id, time_interval]):
        
        # print('Cell {} and period {}'.format(grouping, time_int))
        
        # It may be that there is no burned area in a certain period.
        # In this case, I will print an informative message and I will set the peak month as 0
        if group[var].sum() == 0:
            print('WARNING: no {} in cell {}, time interval {}'.format(var, grouping, time_int))
            m05, m50, m95, length = [np.nan]*4
        else:
            m05, m50, m95, length = calculate_monthDistribution(group, var)
        
        # Appending peak month for this grouping and year to list of dataframes
        var_quant_ts.append(pd.DataFrame({
            grouping_id : grouping,
            time_interval : time_int,
            'q05_month_sECDF': m05,
            'median_month_sECDF': [m50],
            'q95_month_sECDF': m95,
            'q_length_sECDF': length
            }))
        
        
    
    # Converting list of dataframes to single dataframe
    var_quant_ts = pd.concat(var_quant_ts, ignore_index = True)
    
    # Exclude NaNs if asked to
    if exclude_nan == True:
        var_quant_ts = var_quant_ts[~ var_quant_ts['median_month'].isna()].reset_index(drop = True)
        
    return var_quant_ts

def save_fireSeasonCharacteristic(data, data_fp):
    """
    Function that saves the output data to data_fp by previously checking if the
    output directory exists and, if it does not exist, it creates it and then
    saves the file.

    Parameters
    ----------
    data : DataFrame
        Dataframe containing the output data to save.
    data_fp : string
        The filepath where to save the output data.

    Returns
    -------
    None.

    """
    
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
    
    # Total number of periods we are working with
    nperiods = 4
    # Number of years per period
    nyears = 9
    
    # Grid
    grid_id = '30km'
    
    # Directory for input and output files
    dir_fp = '../../data/processed/summary_tables/{}_year_periods/'.format(nyears)

    
    data_fp = dir_fp + 'burnedArea_monthFromPixel_fireCount_assignedCell_exclude3RainyMonths_grid{}_period{}.csv'.format(grid_id, '{}')
    
    # Output filepath
    output_fp = dir_fp + 'fireSeason_measures_exclude3RainyMonths_grid{}_period{}.csv'
    
    
    # --------------- Reading data --------------------------------------------
    
    print('Reading data...')
    data = readFormat_allPeriods(data_fp, nperiods, grid_id)#, thresh = cell_area_thresh)
    
    # ------------------ Processing data ------------------------------------
    
    # Calculating the quantiles per cell and period
    print('Calculating inverse sECDF...')
    inverseCDF = inverse_smoothedECDF_perCell(data, 'id', 'period', 'area', exclude_nan = False)
    
    
    # Subsetting columns of interest
    inverseCDF = inverseCDF[['id', 'period', 'q05_month_sECDF', 'median_month_sECDF', 'q95_month_sECDF', 'q_length_sECDF']]
    
    
    # Freeing memory
    del data
    
    
    # ------------------- Saving data by period -------------------------------
    
    print('Saving output...')    
    
    # Cumulative variable indicating the period we are working on 
    i = 1
    # Looping over the periods, saving the individual files
    for period, data_subset in inverseCDF.groupby('period'):
        
        print('...for period {}'.format(period))
        
        
        # print(data_subset.head(20))
        save_fireSeasonCharacteristic(data_subset, output_fp.format(grid_id, i))
        
        i += 1
    
if __name__ == '__main__':
    main()



