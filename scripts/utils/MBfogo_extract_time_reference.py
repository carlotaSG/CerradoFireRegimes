# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 13:29:18 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Function that returns a time reference (e.g. year or month) reading it from 
the filename.

"""

"""
Imports
"""
import re

"""
Function
"""

def extract_yearMonth(filepath):
    """
    Function that extracts the month and year information from a filename with 
    the format:
        *_<year>-<month>.*
    where year has four digits and month, two.

    Parameters
    ----------
    filepath : str
        Filepath from which to extract the time information.

    Returns
    -------
    year : str
        Four digit year.
    month : str
        Two digit month.

    """
    
    year = re.search(r"_(\d{4})-", filepath).group(1)
    
    month = re.search(r"-(\d{2}).", filepath).group(1)
    
    return year, month

def extract_year(filepath):
    """
    Function that extracts year information from a filename with 
    the format:
        *_<year>.*
    where year has four digits.

    Parameters
    ----------
    filepath : str
        Filepath from which to extract the time information.

    Returns
    -------
    year : str
        Four digit year.
    """
    
    year = re.search(r"_(\d{4}).", filepath).group(1)
    
    return year


