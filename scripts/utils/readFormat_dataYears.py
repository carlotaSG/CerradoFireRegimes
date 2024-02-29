# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 17:04:19 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script originally created for MBfogo_fireSizeCharacteristics.py. It reads the polygon
data and formats it appropriately.

TODO: write proper documentation.
"""
import pandas as pd


import add_gridInfo

def readFormat_dataYears(data_fp, start_year, end_year, grid_id, list_cols,
                         periods = True, nperiods = 4, month = True,
                         fire_polygons = False, min_size = 0, grid_colname = None,
                         drop_geometry = False,
                         thresh = None, incl_cells = None, excl_cells = None):
    
    data = []
    for year in range(start_year, end_year + 1):
        
        # print(year)
        
        # read data for year
        data_year = pd.read_csv(data_fp.format(year))
        
        if fire_polygons == True:
            # Filtering out polygons smaller than x ha
            data_year = data_year[data_year['area'] >= min_size].reset_index(drop = True)
            
            # UPDATE (25/10/2023): with new grids, we do not care about the polygons
            # being on the border of the Cerrado
            # Filtering out border polygons
            # data_year = data_year[data_year['border_pol'] == 0].reset_index(drop = True)
            
            # Renaming cell
            data_year = data_year.rename(columns = {grid_colname : 'id'})
        
        if periods == True:
            if nperiods == 4:
                period_list = [(1985,1993), (1994, 2002), (2003, 2011), (2012, 2020)]
                
                period_labels = ['[1985, 1993]', '[1994, 2002]', '[2003, 2011]', '[2012, 2020]']
                
            # UPDATE 23/08/2023: Included options of just two periods
            elif nperiods == 2:
                period_list = [(1985, 2002), (2003, 2020)]
            
                period_labels = ['[1985, 2002]', '[2003, 2020]']
            
            period_labels = dict(zip(period_labels,
                                     [str(p) for p in period_list]))
            
            period_bins = pd.IntervalIndex.from_tuples(period_list, closed = 'both')
            
            # Assign period
            data_year['period'] = pd.cut(data_year['year'], bins = period_bins).astype(str)
            
            # Replace assigned labels by correct ones
            data_year = data_year.replace({'period':period_labels})
            
            if month == True:
    
                # Selecting columns of interest
                data_year = data_year[['id', 'period', 'year', 'month'] + list_cols]
                
            else:
                # Selecting columns of interest
                data_year = data_year[['id', 'period', 'year'] + list_cols]
                
        else:
            if month == True:
                # Selecting columns of interest
                data_year = data_year[['id', 'year', 'month'] + list_cols]
            else:
                # Selecting columns of interest
                data_year = data_year[['id', 'year'] + list_cols]
        
        
        
        # Appending to data
        data.append(data_year)
        
        # Freeing memory
        del data_year
        
    data = pd.concat(data, ignore_index = True)
    
    # Subsetting those polygons assigned to cells with more than a threshold of
    # their area inside the Cerrado
    data = add_gridInfo.add_gridInfo(data, grid_id, filter_thresh = thresh, excl_cells = excl_cells, incl_cells = incl_cells)
    
    if drop_geometry == True:
        # Dropping geometry column
        data = data.drop(columns = 'geometry')

    return(data)