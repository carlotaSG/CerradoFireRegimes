# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 13:00:10 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that, given a DataFrame or GeoDataFrame, returns the same object with an
extra column containing an ID of the format:
    <year><month>_<dataframe_index>
where the <dataframe_index> is the 0-leftpad id number of the dataframe.
"""

"""
Imports
"""
import math

"""
Function
"""
def create_yearMonthID(dataframe):
    """
    Function that takes a DataFrame of GeoDataFrame,
    assigns a unique identifier to each burned pixel (row) and returns the
    same DataFrame.
    
    Identifier: <year><month>_<dataframe_index>

    Parameters
    ----------
    dataframe : GeoDataFrame
        The read GeoDataFrame.

    Returns
    -------
    dataframe : GeoDataFrame
        The input GeoDataFrame but with an additional column of the burned 
        pixel's unique ID.

    """
    
    # Creating column with the dataframe's index
    dataframe['theindex'] = dataframe.index
    
    #       Number of total digits for the <dataframe_index>
    ndigits = int(math.log10(len(list(dataframe.index)) - 1)) + 1
    
    #       Formatting for the digits to include as many 0 as necessary so all 
    #       of them have the same amount of digits.
    digit_form = "{:0"+str(ndigits)+"d}"
    
    #       Creating column with the unique identifier
    dataframe['ID'] = [str(dataframe.loc[j, 'year']) + 
                       "{:02d}".format(dataframe.loc[j, 'month']) + "_" + 
                       digit_form.format(dataframe.loc[j, 'theindex']) \
                       for j in range(len(dataframe))]
        
    # Cast column as index
    dataframe.index = dataframe['ID']
    
    # Erase temporary column column
    dataframe = dataframe.drop(columns = ['theindex','ID'])
    
    return dataframe

