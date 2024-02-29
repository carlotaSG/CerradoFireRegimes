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
# from statsmodels.distributions.empirical_distribution import ECDF

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
    
    # Calculate month with largest burned area
    peak_month = df.loc[df[var].idxmax(), 'month']
    
    # Calculate how much do we have to add to this value so that this month becomes 6
    offset = 6 - peak_month
    
    # Adding offset to months (relabelling)
    df['month_relabel'] = df['month'] + offset
    
    # Correcting months that now are negative or larger then 12
    df['month_relabel'] = df['month_relabel'].apply(adjustMonth)
    
    return df, offset


def inverseCumPerc(cum_probs, quant):
    
    # Find the index of the inverse ECDF probability closest to quant
    # print(ECDF_probs.sub(quant))
    # print(cum_probs)
    # print(cum_probs.sub(quant))
    # print(cum_probs.sub(quant).abs())
    the_index = cum_probs.sub(quant).abs().idxmin()
    
    # Hence the (relabelled) month with an inverse ECDF closest to quant is
    month_quant = the_index + 1
    
    # if quant == 0.05:
    #     # then we want the closest value to 0.05, but with ECDF greater than 0.05
    #     if ECDF_probs[the_index] < 0.05:
    #         month_quant += 1
        
    #     # Otherwise, we accept the first month assigned
    # elif quant == 0.95:
    #     # In this case, we actually want the month that has ECDF exactly 0.95 or below this value
    #     if ECDF_probs[the_index] > 0.95:
    #         month_quant -= 1
    
    return month_quant

def calculate_empiricalVariance(data, var):
    
    # Calculating the average month in month_relabel
    mean = np.average(data['month_relabel'], weights = data[var])
    
    # Calculate variance
    variance = np.average((data['month_relabel'] - mean)**2, weights = data[var])
    
    # The standard deviation
    std = np.sqrt(variance)
    
    return variance, std


def findQuantile(data, var):
    
    # Order dataframe according to relabelled months
    data = data.sort_values(by = 'month_relabel').reset_index(drop = True)
    
    # # Fit the ECDF to var
    # ecdf = ECDF(data[var])
    
    # # Calculate the inverse ECDF
    # data['ECDF'] = [ecdf(m) for m in range(1, 12+1)]
    
    # Calculate the empirical variance of the distribution and the standard deviation
    mvar, mstd = calculate_empiricalVariance(data, var)
    
    # Calculate the cumulative probability distribution    
    data['cumPerc'] = 100 * (data[var].cumsum() / data[var].sum())
    
    # print(data[['id','month','month_relabel','area','ECDF','cumPerc']])
    
    # Find the quantiles as those (relabelled) months closest to the quantile
    m05 = inverseCumPerc(data['cumPerc'], 5)
    m50 = inverseCumPerc(data['cumPerc'], 50)
    m95 = inverseCumPerc(data['cumPerc'], 95)
    
    # Calculate length of the season. Since in findQuantile we have been conservative,
    # the length of the fire season is the number of months between the two quantiles, 
    # but including the two quantiles
    length = m95 - m05 + 1
    
    # print(m05, m50, m95, length)
    
    return m05, m50, m95, length, mvar, mstd


def calculate_quantilesLength(data, col_var):
    
    # Make sure data is ordered by month
    data = data.sort_values(by = 'month').reset_index(drop = True)
    # print('data')
    # print(data)
    
    # Relabel the months so that peak month is number 6 (distribution is 
    # centered around the peak and in the middle of the fire year)
    data, offset = relabelMonths(data, col_var)
    
    # print(data[['id','month','month_relabel','area']])
    
    
    # Calculate the quantiles obtained when fitting the ECDF, the length of the 
    # season in number of months, and hte variance of the monthly burned area distribution
    m05, m50, m95, length, mvar, mstd = findQuantile(data, col_var)
    
    
    # 4. Convert these quantiles back to the original months
    m05 = adjustMonth(m05 - offset)
    m50 = adjustMonth(m50 - offset)
    m95 = adjustMonth(m95 - offset)
    
    
    return m05, m50, m95, length, mvar, mstd

def inverseCDF_perCell(df, grouping_id, time_interval, var, exclude_nan = False):
    # This function receives a dataframe of the type:
    #   cell/grouping | time_interval | month | variable
    # gropus by cell/grouping and time_interval and calls ranking to find the 
    # legnth of the fire season - the number of months that add up to 80% of the burned area
    # Combinations of cell/grouping and year where variable is - for all months are
    # filled with NA and a message is produced
    
    # List to store the length (number of months) per cell and time_interval in the time series (ts)
    var_quant_ts = []
    
    # Looping over pair grouping - year
    for (grouping, time_int), group in df.groupby( by = [grouping_id, time_interval]):
        
        # print('Cell {}'.format(grouping))
        # print(time_int)
        # It may be that there is no burned area in a certain period.
        # In this case, I will print an informative message and I will set the peak month as 0
        if group[var].sum() == 0:
            print('WARNING: no {} in cell {}, time interval {}'.format(var, grouping, time_int))
            m05, m50, m95, length, mvar, mstd = [np.nan]*6
        else:
            m05, m50, m95, length, mvar, mstd = calculate_quantilesLength(group, var)
        
        # Appending peak month for this grouping and year to list of dataframes
        var_quant_ts.append(pd.DataFrame({
            grouping_id : grouping,
            time_interval : time_int,
            'q05_month': m05,
            'median_month': [m50],
            'q95_month': m95,
            'q_length': length,
            'var_month': mvar,
            'std_month': mstd
            }))
        
        
    
    # Converting list of dataframes to single dataframe
    var_quant_ts = pd.concat(var_quant_ts, ignore_index = True)
    
    # Exclude NaNs if asked to
    if exclude_nan == True:
        var_quant_ts = var_quant_ts[~ var_quant_ts['median_month'].isna()].reset_index(drop = True)
        
    return var_quant_ts

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
    
    # Number of years per period 
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
    
    
    # Minimum area to select cells
    # cell_area_thresh = 75.0
    

    
    
    # --------------- Reading data --------------------------------------------
    
    print('Reading data...')
    data = readFormat_allPeriods(data_fp, nperiods, grid_id)#, thresh = cell_area_thresh)
    
    
    # --------------- Formatting data -----------------------------------------
    
    # # To then identify the peak month
    # # Converting the data from long to wide format
    # data = pd.melt(data, 
    #                id_vars = ['id','%area_crop', 'geometry','period','npolygons','area_T'], 
    #                var_name = 'month', value_name = 'area')
    
    # # Converting month column from string to integer
    # data['month'] = data['month'].apply(int)
    
    
    # --------------- Getting the peak of the fire season ---------------------
    
    print('Calculating the quantiles using the cumulative probability method per cell and period...')
    # data = data[data['period'] == '(2012, 2020)']
    # data = data[data['id'].isin([1528, 1863, 1703, 1851, 1642])]
    inverseCDF = inverseCDF_perCell(data, 'id', 'period', 'area', exclude_nan = False)
    
    
    # Subsetting columns of interest
    inverseCDF = inverseCDF[['id', 'period', 'q05_month', 'median_month', 'q95_month', 'q_length', 'var_month', 'std_month']]
    
    # print(inverseECDF)
    
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
