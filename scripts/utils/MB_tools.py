# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 10:30:25 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script with some variable declaration and other common tasks to perform when
working with MapBiomas data.
"""
import pandas as pd

def dependencies_classLU(collection):
    """
    Three dictionaries with the relations between the sub-categories that are 
    part of a MapBiomas category.

    Returns
    -------
    level_1 : dictionary
        Relations between the most general classes (level_1) and their subcategories
        (level_2).
    level_2 : dictionary
        The sub-sub-categories (level_3) contained in each sub-category (level_2).
    level_3 : dictionary
        The sub-sub-sub-categories (level_4) contained in each 
        sub-sub-category (level_3).

    """
    
    if collection == 5:
    
        level_3 = {
            '19':['39','20','41']
            }
        
        level_2 = {
            '2': ['3','4','5'],
            '18': ['19','36']
            }
        
        level_1 = {
            '1': ['2','9'],
            '10': ['11','12','32','29','13'],
            '14': ['15','18','21'],
            '22': ['23','24','30','25'],
            '26': ['33','31']
            }
    
    elif collection == 6:
        
        level_3 = {
            '19' : ['39','20','40','41'],
            '36' : ['46','47','48']
            }
        
        level_2 = {
            '18': ['19','36']
            }
        
        level_1 = {
            '1'  : ['3','4','5','49'],
            '10' : ['11','12','32','29','13'],
            '14' : ['15','18','9','21'],
            '22' : ['23','24','30','25'],
            '26' : ['33','31']
            }
        
    elif collection == 8:
        
        level_3 = {
            '19' : ['39','20','40','41','62'],
            '36' : ['35','46','47','48']
            }
        
        level_2 = {
            '18': ['19','36']
            }
        
        level_1 = {
            '1'  : ['3','4','5','6','49'],
            '10' : ['11','12','32','29','13','50'],
            '14' : ['15','18','9','21'],
            '22' : ['23','24','30','25'],
            '26' : ['33','31']
            }
    
    return level_1, level_2, level_3

def standard_LU(collection):
    """
    List of classes that I consider standard as strings: the main natural classes
    and the different types of agricultural classes

    Returns
    -------
    standard_landuses : list of strings
        List of land-use classes.

    """
    if collection == 5:
        standard_landuses = [3, 4, 9, 12, 15, 39, 20, 41, 36, 21]
    elif collection == 6:
        standard_landuses = [3, 4, 11, 12, 15, 18, 21, 22]
    
    # Do I need this in strign mode?
    standard_landuses = [str(i) for i in standard_landuses]
    
    return standard_landuses

def main_LU_c5():
    """
    List of classes that I consider principal as strings: the main natural classes
    and the main agricultural classes

    Returns
    -------
    main_landuses : list of strings
        List of land-use classes.

    """
    main_landuses = [3, 4, 12, 18, 15, 21, 22]
    
    main_landuses = [str(i) for i in main_landuses]
    
    return main_landuses

def main_cutoff_LU():
    """
    List of classes that I consider principal, except some sub-classes of land-use
    as strings: the main natural classes and the main agricultural classes

    Returns
    -------
    main_cutoff_LU : list of strings
        List of land-use classes.

    """
    main_landuses = [3, 4, 12, 18, 15]
    
    main_landuses = [str(i) for i in main_landuses]
    
    return main_landuses

def select_LU(list_columns, list_LU, df):
    """
    Given a list of land-use classes and columns to keep, and a dataframe in wide
    format (land-uses are the columns of the dataframe), it selects only those 
    columns in list_LU    

    Parameters
    ----------
    list_columns : list of strings
        List of non land-use columns to keep.
    list_LU : list of strings
        List of land-use columns to keep.
    df : DataFrame
        DataFrame from which to keep columns.

    Returns
    -------
    df : DataFrame
        The subset DataFrame.

    """
    #Columns to keep
    keep_columns = list_columns + list_LU
    
    df = df[keep_columns]
    
    return df

def wide_to_long_LU(df, columns_index, val_name):
    """
    Receives a dataframe in wide format and casts it to long format:
    with columns_index to remain and the land-use columns to be cast into 
    long format:
    columns_index | category_LU | val_name    

    Parameters
    ----------
    df : DataFrame
        DataFrame to cast into long format.
    columns_index : List of strings
        List of columns to keep in wide format.
    val_name : string
        Name to be given to the column that will contain the content of the 
        land-use columns.

    Returns
    -------
    df : DataFrame
        DataFrame in long format.

    """
    df = pd.melt(df, id_vars = columns_index, var_name= 'category_LU', value_name=val_name)
    
    return df

def pixel_to_area_km2(df, id_columns):
    """
    Function that multiplies all elements in dataframe by a constant, except
    the list of id_columns: converts all values that are number of pixels to
    area in km2.

    Parameters
    ----------
    df : DataFrame
        DataFrame whose values are to be transformed from number of pixels to
        area (Landsat spatial resolution.
    id_columns : list of strings
        List of columns that do not have to be covnerted to area.

    Returns
    -------
    df : DataFrame
        DataFrame with values transformed to area.

    """   
    area_conversion = 30 * 30 * 1e-06
    df.loc[:, ~df.columns.isin(id_columns)] = df.loc[:, ~df.columns.isin(id_columns)] * area_conversion
    
    return df

def fill_classes(df, collection):
    """
    Function that calculates the values for those MapBiomas LU classes 
    that have subclasses and hence do not appear directly in the raster images.    

    Parameters
    ----------
    df : DataFrame
        DataFrame with some empty categories and sub-categories.

    Returns
    -------
    df : DataFrame
        DataFrame with no empty categories and sub-categories.

    """
    # Getting dependencies relations
    level_1, level_2, level_3 = dependencies_classLU(collection)
    
    # LEVEL 3
    for lu_class in level_3:
        df[lu_class] = df[level_3.get(lu_class)].sum(axis=1)
        
    # LEVEL 2
    for lu_class in level_2:
        df[lu_class] = df[level_2.get(lu_class)].sum(axis=1)
        
    # LEVEL 3
    for lu_class in level_1:
        df[lu_class] = df[level_1.get(lu_class)].sum(axis=1)
    
    return df

def monthCode():
    monthCode = {
        'Jan' : '01',
        'Feb' : '02',
        'Mar' : '03',
        'Apr' : '04',
        'May' : '05',
        'Jun' : '06',
        'Jul' : '07',
        'Aug' : '08',
        'Sep' : '09',
        'Oct' : '10',
        'Nov' : '11',
        'Dec' : '12'
        }
    
    return monthCode