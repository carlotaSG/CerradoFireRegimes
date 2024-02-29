# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 11:45:43 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculates the fire season using:
    - Definition 1: define a fire year as a 12 year period where month number 
    is the month with more fires in a calendar year. All other months are 
    relabelled accordingly, with Januray going after December. Then, the fire
    season is delimited as the period containing a certain interquantile.
    
This script is called as a function from other script. The function called is
calculate_season.

"""

"""
Imports
"""
import pandas as pd

"""
Functions
"""

def calculate_transformation(data):
    """
    Function that applieds Definition 1 as fire season definition.
    Gets the list of cells along with the month when most fires occurred, it 
    subtracts this month value from 6 to assign it as month 0, the difference 
    is the amount of months the fire polygon's months must be transformed by.

    Parameters
    ----------
    data : DataFrame
        Containing the list of cells and their month with more fires.

    Returns
    -------
    data : DataFrame
        Returns input dataframe along with a column 'trasnformation' indicating
        the number by which fire polygon's months must be transformed.

    """

    # Calculating translation
    data['transformation'] = 6 - data['month']
    
    
    return data

def adjustMonth(month, transformation):
    """
    Once the initial set of months are transformed by adding the transformation 
    value, seom months may be negative and others may be larger than 12. Hence, 
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
    
def apply_translation(data, month_info):
    """
    Function that, given the information on the month that has to be used as a 
    reference and the list of months, transforms all months according to the
    fire season definition and the month_info reference.

    Parameters
    ----------
    data : DataFrame
        The DataFrame containing the list of polygons assigned to the cell and
        the month to which they belong.
    month_info : DataFrame
        Contains the list of cells with the corresponding month that is to be 
        used as a reference.

    Returns
    -------
    new_data : DataFrame
        Contains the same information as the variable data but with the month 
        transformed according to the fire year definition used and the month reference.

    """
    
    # First, calculate month transformation
    trans_data = calculate_transformation(month_info)
    
    # Grouping the data by 'id' and applying corresponding transformation of month 
    # column
    new_data = []
    for unit, data_group in data.groupby(by = 'id'):
        
        # Getting the transformation for this unit 
        trans_group = trans_data.loc[trans_data['id'] == unit, 'transformation'].item()
        
        data_group['month_trans'] = data_group['month'] + trans_group

        # Adjusting month_trans label to be within 1 and 12
        data_group['month_trans'] = data_group['month_trans'].apply(adjustMonth, transformation = trans_group)
        
        
        # Appending modified group to list of dataframes
        new_data.append(data_group)
        
    
    # Converting list of lists to unique dataframe
    new_data = pd.concat(new_data, ignore_index = True)
    
    return new_data

def calculate_season(data, quant, month_info):
    """
    Function that calculates the fire season. The start, end and length of the
    season.

    Parameters
    ----------
    data : DataFrame
        The DataFrame containing the list of polygons assigned to the cell and
        the month to which they belong.
    quant : float
        The percentage to be used as interquantile.
    month_info : DataFrame
        Contains the list of cells with the corresponding month is to be 
        used as a reference.

    Returns
    -------
    season : DataFrame
        Contains the list of cells with their corresponding start, end and 
        length of the fire season.

    """
    
    # TODO: add definition option
    
    # First, transform the data according to definition
    transformed_data = apply_translation(data, month_info)
    
    # Start of fire season in transformed label
    start_season = transformed_data[['id', 'month_trans']].groupby(by = 'id').quantile(quant[0], interpolation = 'nearest').reset_index()
    
    # End of fire season in transformed label
    end_season = transformed_data[['id', 'month_trans']].groupby(by = 'id').quantile(quant[1], interpolation = 'nearest').reset_index()
    
    # Merging datasets
    season = pd.merge(start_season, end_season, how = 'outer', on = 'id',
                      suffixes = ('_start', '_end'))
    
    # Calculating season length
    season['length'] = season['month_trans_end'] - season['month_trans_start'] + 1
    
    # Merging to month_info
    season = pd.merge(season, month_info[['id', 'transformation']], on = 'id')
    
    # Getting the start and end of the fires season in 'original month labels'
    season['month_start'] = season['month_trans_start'] - season['transformation']
    season['month_end'] = season['month_trans_end'] - season['transformation']
    
    season['month_start'] = [adjustMonth(season.loc[i, 'month_start'], season.loc[i, 'transformation']) for i in range(len(season))]
    season['month_end'] = [adjustMonth(season.loc[i, 'month_end'], season.loc[i, 'transformation']) for i in range(len(season))]
    
    # Selecting relevant columns
    season = season[['id', 'month_start', 'month_end', 'length']]
    
    del transformed_data, start_season, end_season
    
    return season
        