# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 17:18:35 2022

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that finds the peak of the fire season for a cell in a period of time 
(provided by the user as time_interval). It also contains a function that 
simply finds the month with larges var value.


"""

"""
Imports
"""
import pandas as pd
import numpy as np

"""
Functions
"""
def find_peak(df, var):
    # Function that receives a df of the form
    #   month | var
    # and finds the month for which var is maximum
    
    peak = df.loc[df[var].idxmax(), 'month']
    
    return peak
    

def find_peak_perCell(df, grouping_id, time_interval, var, exclude_nan = False):
    # This function receives a dataframe of the type:
    #   cell/grouping | time_interval | month | variable
    # gropus by cell/grouping and time_interval and calls find_peak to find the month
    # peak of variable per year and cell/grouping
    # Combinations of cell/grouping and year where variable is - for all months are
    # filled with NA and a message is produced
    
    # List to store the peak month per cell and time_interval in the time series (ts)
    var_peak_ts = []
    
    # Looping over pair grouping - year
    for (grouping, time_int), group in df.groupby( by = [grouping_id, time_interval]):
        # It may be that there is no burned area in a certain year.
        # In this case, I will print an informative message and I will set the peak month as 0
        if group[var].sum() == 0:
            print('WARNING: no {} in cell {}, time interval {}'.format(var, grouping, time_int))
            peak_month = np.nan
        else:
            peak_month = find_peak(group, var)
        
        # Appending peak month for this grouping and year to list of dataframes
        var_peak_ts.append(pd.DataFrame({
            grouping_id : grouping,
            time_interval : time_int,
            'peak_month' : [peak_month]
            }))
    
    # Converting list of dataframes to single dataframe
    var_peak_ts = pd.concat(var_peak_ts, ignore_index = True)
    
    # Exclude NaNs is asked to
    if exclude_nan == True:
        var_peak_ts = var_peak_ts[~ var_peak_ts['peak_month'].isna()].reset_index(drop = True)
        
    return var_peak_ts
            
            
        
    
    

"""
MAIN
"""

