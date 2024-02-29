# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 13:41:29 2021

@author: scat8298

Function that matches filenames based on some criteria like month and year
from two different sets of filename lists and returns
a list of tuple, where each tuple is (filename_list1, filename_list2).
"""

"""
Imports
"""
import pandas as pd
from functools import reduce
from MBfogo_extract_time_reference import extract_yearMonth, extract_year

import sys

"""
Functions
"""
def createDF_yearMonth(list_fp):
    # Function that turns a list of filpeaths into a dataframe of
    # columns:
        # fn | yearMonth
    
    # Cast first list into a dataframe 
    list_fp = pd.DataFrame({'fn' : list_fp})
    
    # Add year and months as columns
    yearMonth = list(map(extract_yearMonth, list_fp['fn']))
    
    list_fp['yearMonth'] = [element[0] + element[1] for element in yearMonth]
    #list_fp['month'] = [element[1] for element in yearMonth]

    return list_fp

def createDF_year(list_fp):
    # Function that turns a list of filpeaths into a dataframe of
    # columns:
        # fn | year
    
    # Cast first list into a dataframe 
    list_fp = pd.DataFrame({'fn' : list_fp})
    
    # Add year and months as columns
    year = list(map(extract_year, list_fp['fn']))
    
    list_fp['year'] = year # [element[0] + element[1] for element in year]
    #list_fp['month'] = [element[1] for element in yearMonth]

    return list_fp


def matchingFileNames_monthYear(list1, list2):
    # Functiont that gets two lists and returns a list of tuples
    # with all the elements in list1, and the corresponding elements with same year and month
    # in list2 in second position in tuple 
    
    df1 = createDF_yearMonth(list1)
    
    df2 = createDF_yearMonth(list2)
    
    df = df1.merge(df2, how='left', on='yearMonth')
    
    list_tuples = list( zip( df['fn_x'], df['fn_y']))
    
    return list_tuples

def matchingFileNames_year(list_of_fp_list):#list1, list2):
    # Functiont that gets two lists and returns a list of tuples
    # with all the elements in list1, and the corresponding elements with same year 
    # in list2 in second position in tuple 
    
    list_df = []
    for list_fp in list_of_fp_list:
        
        # Create dataframe
        df = createDF_year(list_fp)
        
        # Sort dataframe by year
        df = df.sort_values('year')
        
        # Append dataframe to list of dataframes
        list_df.append(df)
        
    # Checking if list of years is the same for all dataframes
    list_years = [df['year'].tolist() for df in list_df]
    result = all(element == list_years[0] for element in list_years)
    if result == False:
        print('WARNING: matchingFileNames.matchingFileNames_year')
        print('The lists of filepaths passed to function do not cover the same set of years.')
        print('Exiting program.')
        sys.exit()
    
    #list_fp = [df['fn'].tolist() for df in list_df]
    
    # Merge dataframes
    df_merged = reduce(lambda  left, right: pd.merge(left, right, on=['year'],
                                            how='outer'), list_df)
    
    #print(df_merged.columns)
    
    # Drop the 'year' column as it is not needed anymore
    df_merged = df_merged.drop('year', axis=1)
    
    # Create list of tuples with same year filepaths
    list_tuples = df_merged[df_merged.columns].apply(tuple, axis=1).tolist()
    
    # print(output)
    
    # list_tuples = list(zip([element ]))
    
    # print(list_tuples)
    # df1 = createDF_year(list1)
    
    # df2 = createDF_year(list2)
    
    # df = df1.merge(df2, how='left', on='year')
    
    # list_tuples = list( zip( df['fn_x'], df['fn_y']))
    
    return list_tuples
    
    
