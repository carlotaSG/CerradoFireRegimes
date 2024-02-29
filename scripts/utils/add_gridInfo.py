# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 16:17:15 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that adds the geometric information of a grid to a set of data. Apart from
the geoemtry, it returns the %area_crop, the percentage of a grid cell inside the
Cerrado. It also allows the user to pass a %area_crop threshold, in which case
the script returns only those cells with an area equal or greater than the threshold.


"""

"""
Imports
"""
import geopandas as gpd

"""
Functions
"""
def add_gridInfo(data, grid_id, filter_thresh = None, excl_cells = None, incl_cells = None):
    # Function that merges the geometry information to a grid data and, optionally, selects only those cells with %area_crp above a certain threshld
    
    # Reading grid
    grid = gpd.read_file('../../data/shapes/cerrado_grid_{}_cropped.shp'.format(grid_id))
    grid = grid[['id','%area_crop','geometry']]
    
    # Merging geographic information with data
    data = grid.merge(data, on = 'id', how = 'right')
    
    if (filter_thresh != None) & (excl_cells != None) & (incl_cells != None):
        # Include only cells with area greater than threshold
        data = data[((data['%area_crop'] >= filter_thresh) | data['id'].isin(incl_cells)) & ~data['id'].isin(excl_cells)]
        
    elif (filter_thresh != None) & (excl_cells != None) & (incl_cells != None):
        # Include only cells with area greater than threshold
        data = data[((data['%area_crop'] >= filter_thresh) | data['id'].isin(incl_cells)) & ~data['id'].isin(excl_cells)]
    elif (filter_thresh != None) & (excl_cells == None) & (incl_cells == None):
        # Include only cells with area greater than threshold
        data = data[data['%area_crop'] >= filter_thresh]
    
    return(data)