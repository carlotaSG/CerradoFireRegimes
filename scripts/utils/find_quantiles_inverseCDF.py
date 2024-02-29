# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 11:50:13 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk
"""
"""
Imports
"""
import numpy as np
import pandas as pd

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


def findQuantile(data, var):
    
    # Order dataframe according to relabelled months
    data = data.sort_values(by = 'month_relabel').reset_index(drop = True)
    
    # # Fit the ECDF to var
    # ecdf = ECDF(data[var])
    
    # # Calculate the inverse ECDF
    # data['ECDF'] = [ecdf(m) for m in range(1, 12+1)]
    
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
    
    return m05, m50, m95, length


def calculate_quantilesLength(data, col_var):
    
    # Make sure data is ordered by month
    data = data.sort_values(by = 'month').reset_index(drop = True)
    
    # Relabel the months so that peak month is number 6 (distribution is 
    # centered around the peak and in the middle of the fire year)
    data, offset = relabelMonths(data, col_var)
    
    # print(data[['id','month','month_relabel','area']])
    
    # Calculate the quantiles obtained when fitting the ECDF
    m05, m50, m95, length = findQuantile(data, col_var)
    
    
    # 4. Convert these quantiles back to the original months
    m05 = adjustMonth(m05 - offset)
    m50 = adjustMonth(m50 - offset)
    m95 = adjustMonth(m95 - offset)
    
    
    return m05, m50, m95, length

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
        # It may be that there is no burned area in a certain period.
        # In this case, I will print an informative message and I will set the peak month as 0
        if group[var].sum() == 0:
            print('WARNING: no {} in cell {}, time interval {}'.format(var, grouping, time_int))
            m05, m50, m95, length = [np.nan]*4
        else:
            m05, m50, m95, length = calculate_quantilesLength(group, var)
        
        # Appending peak month for this grouping and year to list of dataframes
        var_quant_ts.append(pd.DataFrame({
            grouping_id : grouping,
            time_interval : time_int,
            'q05_month': m05,
            'median_month': [m50],
            'q95_month': m95,
            'q_length': length
            }))
        
        
    
    # Converting list of dataframes to single dataframe
    var_quant_ts = pd.concat(var_quant_ts, ignore_index = True)
    
    # Exclude NaNs if asked to
    if exclude_nan == True:
        var_quant_ts = var_quant_ts[~ var_quant_ts['median_month'].isna()].reset_index(drop = True)
        
    return var_quant_ts