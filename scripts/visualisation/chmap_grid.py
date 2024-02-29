# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 16:40:26 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that 

Update 02/02/2022: Introduced modification to only plot cells with more than 70%
                   of their area with in the Cerrado
"""

"""
Imports
"""
import pandas as pd
import geopandas as gpd
import numpy as np


import sys
sys.path.append('../../scripts/')

import chmap
import chmap_byMonth
import hist
import tools
import MB_legend_codes as mblegend

"""
Functions
"""
def chmap_burnedArea(grid_size, temp):
    
    # Relevant data files
    # Data file
    data_fp = '../data/processed/grid_{}km/MBfogo_burnedArea_grid{}km_cerrado_2019.csv'.format(grid_size, grid_size)
    # Grid SHP file
    grid_fp = '../data/grids/cerrado_grid_{}km_clipped.shp'.format(grid_size)
    
    # Reading data file
    data = pd.read_csv(data_fp)
    
    # Update 02/02/2022: Filtering out cells with less than 70% of their area
    # within the Cerrado
    data = data[data['%area_crop'] >= 70.0]
    
    month_columns = [str(i) for i in range(1, 12 + 1)]
        
    # Depending on temp, summarise data and transform number of pixels to area
    if temp == 'm':
        # If resolution is monthly, we do not summarise the data, only convert
        # to area in km2
        data.loc[:, month_columns] = data.loc[:, month_columns] * 30 * 30 * 1e-06
        
        # Getting minimum and maximum burned area values, to be used in monthly plot
        min_val = data[month_columns].min()
        min_val = min_val.min()
        
        max_val = data[month_columns].max()
        max_val = max_val.max()
        
    elif temp == 'y':
        # Summarise
        data['total'] = data[month_columns].sum(axis = 1)
        
        # Select only relevant columns
        data = data[['id', 'total']]
        
        # Convert number of pixels to area
        data['total'] = data['total'] * 30 * 30 * 1e-06
        
    # Read grid
    grid = gpd.read_file(grid_fp)    
        
    # Merging geometry and data files by id column
    grid = grid.merge(data, on = 'id', how = 'right')
    
    
    # Plot figure depending on temp
    if temp == 'y':
        # Choropleth map of total burned area
        chmap.chmap(grid, 'total', 'Burned area (km2)', 
                    title = 'Total burned area in the Cerrado in 2019. Grid of {} km'.format(grid_size),
                    filename = 'MBfogo_burnedArea_map_grid{}km.png'.format(grid_size))
        
        # Calculating percentage of burned area for each cell - counting only cell area inside the Cerrado
        grid['percArea'] = 100 * grid['total'] / (grid['area']*1e-6)
        
        # Choropleth map of of burned area percentage
        chmap.chmap(grid, 'percArea', 'Burned area (km2)', 
                    title = 'Percentage burned area in the Cerrado in 2019. Grid of {} km'.format(grid_size),
                    filename = 'MBfogo_burnedAreaPerc_map_grid{}km.png'.format(grid_size))

        # # Filtering out cells with no burned area
        # # grid = grid[grid['percArea'] > 0]
        # print(len(grid))
        # # Histogram
        # hist.hist(grid, 'percArea', 
        #           title = 'Histrogram of percentage of burned area per grid cell (Cerrado, 2019). \n Grid of {}km'.format(grid_size),
        #           bins = 20,
        #           xlabel = 'Burned area percentage', ylabel = 'Number of cells',
        #           filename = 'MBfogo_burnedArea_hist_percArea_grid{}km.png'.format(grid_size))
        
        # # Custom bins
        # bin_edges = [0.0, 0.0, 0.1, 1.0, 10.0, 25.0, 50., 100., 200., 400., 600., 1250.]
        
        hist.hist(grid, 'total', 
                  title = 'Histogram of burned area per grid cell 2019 \n Grid of {}km'.format(grid_size),
                  bins = 20,
                  xlabel = 'Burned area', ylabel = 'Number of cells',
                  filename = 'MBfogo_burnedArea_hist_grid{}km.png'.format(grid_size))
    elif temp == 'm':
        chmap_byMonth.chmap_byMonth(grid, 'Burned area (km2)',
                                    title = 'Burned area per month in the Cerrado in 2019. \n Grid of {} km'.format(grid_size),
                                    filename = 'MBfogo_burnedArea_monthlyMap_grid{}km.png'.format(grid_size),
                                    min_val = min_val, max_val = max_val)
    
    
def chmap_numPolygons(grid_size, temp, criteria):

    # Relevant data files
    # Data file
    if temp == 'y':
        # data_fp = '../data/processed/MBfogo_polygon_grid{}km_counts_cerrado_2019.csv'.format(grid_size)
        # Update 02/02/2022: Read the file that now contains this information 
        # (which is also updated to not include small polygons nor border polygons,
        # nore cells with less than 70% area inside Cerrado)
        data_fp = '../data/processed/grid_{}km/summary_tables/summary_fire-traits_grid50km_2019.csv'.format(grid_size)
        
        # Reading data file
        data = pd.read_csv(data_fp, usecols = ['id','npolygons'])
        
    elif temp == 'm':
        data_fp = '../data/processed/MBfogo_polygon_grid{}km_monthly-counts_cerrado_2019.csv'.format(grid_size)
        
        # Reading data file
        data = pd.read_csv(data_fp)
        
    # Grid SHP file
    grid_fp = '../data/grids/cerrado_grid_{}km_cropped.shp'.format(grid_size)    
    
    
    
    # Read grid
    grid = gpd.read_file(grid_fp)    
        
    # print(grid)
    
    
    month_columns = [str(i) for i in range(1, 12 + 1)]
    # month_columns = list(range(1, 12 + 1))
        
    # Update 02/02/2022   
    # if criteria == 1:
    #     variable = 'npol_inters'
    # elif criteria == 2:
    #     variable = 'npol_assign'
        
    # Depending on temp, summarise data and merge to grid
    # if temp == 'y':
    #     data = data[['cell_ID', variable]]
    if temp == 'm':
        # Converting month column from integer to string, as is format required in chmap_byMonth
        data['month'] = [str(element) for element in data['month']]
        # In this case we transform the long format to wide format (one column per month)
        data = data[['cell_ID', 'month', 'npolygons']].pivot_table(index = ['cell_ID'],
                                                                columns = 'month',
                                                                values = 'npolygons').reset_index()
        
        # print(data)
        # print(data.columns)
        # Extract minimum value and maximum value. Doing this so that all monthly maps
        # share the same legend
        min_val = data[month_columns].min()
        min_val = min_val.min()
        
        max_val = data[month_columns].max()
        max_val = max_val.max()
        
        
        
        
        
        
    # Merging geometry and data files by id column
    grid = grid.merge(data, on = 'id', how = 'right')
    
    grid = grid.fillna(0)
    
    
    # print(grid)
    
    
    # Plot figure depending on temp
    if temp == 'y':
        chmap.chmap(grid, 'npolygons', 'Number of fire polygons', 
                    title = 'Number of fire polygons per cell in the Cerrado in 2019. \n Grid of {} km'.format(grid_size),
                    filename = 'MBfogo_numPolygons_map_grid{}km_criteria{}.png'.format(grid_size, criteria))
        
        hist.hist(grid, 'npolygons', 
            title = 'Histogram of number of fire polygons per grid cell (Cerrado, 2019) \n Grid of {}km'.format(grid_size),
            bins = 30,
            rangeb = (1, grid['npolygons'].max()),
            xlabel = 'Number of fire polygons', ylabel = 'Number of cells',
            filename = 'MBfogo_numPolygons_hist_grid{}km.png'.format(grid_size))
        
    elif temp == 'm':
        chmap_byMonth.chmap_byMonth(grid, 'Number of fire polygons', 
                                    title = 'Number of fire polygons per month and cell in the Cerrado in 2019/ \n Grid of {} km and criteria {}'.format(grid_size, criteria),
                                    filename = 'MBfogo_numPolygons_monthlyMap_grid{}km_critera{}.png'.format(grid_size, criteria),
                                    min_val = min_val, max_val = max_val)
        
    # Now let's plot the same graph but only distinguishing if cell is empty (1) or not (0)
    
    # if temp == 'y':
    #     # Add one column with 0 and 1
    #     grid['nofires'] = np.where(grid[variable] == 0.0, 1, 0)
    #     print(grid)
    #     chmap.chmap(grid, 'nofires', 'No fires (1), at least one fire (0)', 
    #                 title = 'Cells with no fires in the Cerrado in 2019. \n Grid of {} km'.format(grid_size),
    #                 filename = 'MBfogo_noFires_map_grid{}km_criteria{}.png'.format(grid_size, criteria),
    #                 categorical = True)

def chmap_sizeQuantile(grid_size, temp, criteria):
    
    # First select data file
    # Update 02/02/2022: this data file is obsolete, isntead we use
    # data_fp = '../data/processed/MBfogo_polygon_grid{}km_cerrado_2019_criteria0{}.csv'.format(grid_size, criteria)
    data_fp = '../data/processed/MBfogo_c10_cerrado_polygons/MBfogo_c10_cerrado_polygons_2019.csv'
    
    # Then, grid file
    grid_fp = '../data/grids/cerrado_grid_{}km_cropped.shp'.format(grid_size)
    
    # Read both files
    data = pd.read_csv(data_fp)
    # Update 02/02/2022 filter out small polygons and border polygons
    data = data[(data['area'] > 0.009) & (data['border_pol'] == 0)]
    
    
    grid = gpd.read_file(grid_fp)
    # Update 02/02/2022: Select only grid cells with >= 70% area inside Cerrado
    grid = grid[grid['%area_crop'] >= 70.0]
    
    month_columns = [str(i) for i in range(1, 12 + 1)]

    
    # Set cell_ID column name depending on critera
    # First, get name based on criteria
    # if criteria == 1:
    #     variable = 'cell_inters'
    # elif criteria == 2:
    #     variable = 'cell_assign'
        
    # # Then, change to 'cell_ID'
    # data = data.rename(columns = {variable : 'cell_ID'})
    # print(data)
    
    
    # # Pre-arragning data if temporal interval is 'month'
    # if temp == 'm':
    #     # Converting month column from integer to string, as is format required in chmap_byMonth
    #     data['month'] = [str(element) for element in data['month']]
    #     # In this case we transform the long format to wide format (one column per month)
    #     data = data[['cell_ID', 'month', 'area']].pivot_table(index = ['cell_ID'],
    #                                                             columns = 'month',
    #                                                             values = 'area').reset_index()
        
        
    # print(data)


    
    # List of quantiles for which we want to calculate the corresponding areas
    quantiles = [0.5, 0.95, 1.0]
    
    # Generate one plot per quantile
    for q in quantiles:
        
    
    
        # Depending on the temporal interval, calculate quantiles by month or for
        # the whole year
        if temp == 'y':
            data_quant = data[['cellID_50km', 'area']].groupby(by = 'cellID_50km').quantile(q).reset_index()
            # print(data_quant)
            
        elif temp == 'm':
            # First, calculate quantile for each month-cell_ID combination
            data_quant = data.groupby(by = ['month', 'cellID_50km']).quantile(q).reset_index()
            
            # Then, pivot to wide-format
            # Converting month column from integer to string, as is format required in chmap_byMonth
            data_quant['month'] = [str(element) for element in data_quant['month']]
            # In this case we transform the long format to wide format (one column per month)
            data_quant = data_quant[['cellID_50km', 'month', 'area']].pivot_table(index = ['cellID_50km'],
                                                                    columns = 'month',
                                                                    values = 'area').reset_index()
            
            # print(data_quant)
            min_val = data_quant[month_columns].min()
            min_val = min_val.min()
            
            max_val = data_quant[month_columns].max()
            max_val = max_val.max()
            
        # Merge to grid
        grid_quant = grid.merge(data_quant, left_on = 'id', right_on = 'cellID_50km', how = 'left').copy()
        
        
        # PLot figure depending on temp
        if temp == 'y':
            chmap.chmap(grid_quant, 'area', 'Area (km2)', 
                        title = 'Fire polygon area for quantile {}% for the Cerrado in 2019. \n Grid of {} km'.format(int(100*q), grid_size),
                        filename = 'MBfogo_polArea_quantile{}_map_grid{}km.png'.format(int(100*q), grid_size))
        elif temp == 'm':
            chmap_byMonth.chmap_byMonth(grid_quant, 'Area (km2)', 
                                        title = 'Fire polygon area for quantile {}% for the Cerrado in 2019. \n Grid of {} km and criteria {}'.format(int(100*q), grid_size, criteria),
                                        filename = 'MBfogo_polArea_quantile{}_monthlyMap_grid{}km_criteria{}.png'.format(int(100*q), grid_size, criteria),
                                        min_val = min_val, max_val = max_val)
def select_monthMoreFires_perCell(data):
    
    data = data[['month', 'cell_ID']].groupby(by = ['cell_ID','month']).size().reset_index(name = 'npol')
    
    # Subsetting the dataset by selecting the month that has the largest number of fires per cell
    data_subset = data.iloc[data[['cell_ID','npol']].groupby( by = ['cell_ID']).idxmax().reset_index()['npol'], :]
    
    
    return data_subset
    
def chmap_seasonality_moreFires(grid_size, temp, criteria):    
    
    data_fp = '../data/processed/MBfogo_polygon_grid{}km_cerrado_2019_criteria0{}.csv'.format(grid_size, criteria)
    
    grid_fp = '../data/grids/cerrado_grid_{}km_cropped.shp'.format(grid_size)

    # Reading files
    data = pd.read_csv(data_fp)
    grid = gpd.read_file(grid_fp)
    
        # Set cell_ID column name depending on critera
    # First, get name based on criteria
    if criteria == 1:
        variable = 'cell_inters'
    elif criteria == 2:
        variable = 'cell_assign'
        
    # Then, change to 'cell_ID'
    data = data.rename(columns = {variable : 'cell_ID'})
    
    data_subset = select_monthMoreFires_perCell(data)
    
    # data = data[['month', 'cell_ID']].groupby(by = ['cell_ID','month']).size().reset_index(name = 'npol')
    
    # # Subsetting the dataset by selecting the month that has the largest number of fires per cell
    # data_subset = data.iloc[data[['cell_ID','npol']].groupby( by = ['cell_ID']).idxmax().reset_index()['npol'], :]
    
    grid = grid.merge(data_subset, left_on = 'id', right_on = 'cell_ID', how = 'outer')
    
    chmap.chmap(grid, 'month', 'Month', 
                title = 'Month with more fires in the Cerrado in 2019. \n Grid of {} km and criteria {}'.format(grid_size, criteria),
                filename = 'MBfogo_moreFiresMonth_map_grid{}km_criteria{}.png'.format(grid_size, criteria),
                categorical = True, colmap = 'viridis')
    
    hist.hist(grid, 'month', 
            title = 'Histogram of number of cells against month with more fires (peak of fire season) (Cerrado, 2019) \n Grid of {}km'.format(grid_size),
            bins = np.arange(1, 13+1) -0.5,
            # rangeb = (1, 12),
            xlabel = 'Month', ylabel = 'Number of cells',
            filename = 'MBfogo_moreFiresMonth_hist_grid{}km.png'.format(grid_size))

def chmap_seasonality_largestFire(grid_size, temp, criteria):    
    
    data_fp = '../data/processed/MBfogo_polygon_grid{}km_cerrado_2019_criteria0{}.csv'.format(grid_size, criteria)
    
    grid_fp = '../data/grids/cerrado_grid_{}km_cropped.shp'.format(grid_size)

    # Reading files
    data = pd.read_csv(data_fp)
    grid = gpd.read_file(grid_fp)
    
    # Set cell_ID column name depending on critera
    # First, get name based on criteria
    if criteria == 1:
        variable = 'cell_inters'
    elif criteria == 2:
        variable = 'cell_assign'
        
    # Then, change to 'cell_ID'
    data = data.rename(columns = {variable : 'cell_ID'})
    
    # print(data)
    # Subsetting the dataset by selecting for each cell, the month with the largest fire
    data_subset = data.iloc[data[['cell_ID','area']].groupby(by = ['cell_ID']).idxmax().reset_index()['area'], :]
    data_subset = data_subset[['cell_ID', 'month']].reset_index(drop = True)
    print(data_subset)

    grid = grid.merge(data_subset, left_on = 'id', right_on = 'cell_ID', how = 'outer')
    
    chmap.chmap(grid, 'month', 'Month', 
                title = 'Month with largest fire per cell in the Cerrado in 2019. \n Grid of {} km and criteria {}'.format(grid_size, criteria),
                filename = 'MBfogo_largestFireMonth_map_grid{}km_criteria{}.png'.format(grid_size, criteria),
                categorical = True, colmap = 'viridis')
    
    hist.hist(grid, 'month', 
            title = 'Histogram of number of cells against month with largest fire (Cerrado, 2019) \n Grid of {}km'.format(grid_size),
            bins = np.arange(1, 13+1) -0.5,
            # rangeb = (1, 12),
            xlabel = 'Month', ylabel = 'Number of cells',
            filename = 'MBfogo_largestFireMonth_hist_grid{}km.png'.format(grid_size))

def season_code(start, end):
    
    # First, get dictionary of month codes
    month_codes = tools.get_monthName()
    
    # Convert dictionary values to list
    month_codes = list(month_codes.values())
    
    # Convert both the start and end series of floats to integers and then to strings
    start = start.apply(int)
    end = end.apply(int)
    
    season_code = []
    for i in range(len(start)):
        
        # If range of month labels is ascending
        if end[i] >= start[i]:
            # Select the range of codes
            season_code.append(''.join(month_codes[start[i] - 1 : end[i]]))
        elif end[i] < start[i]:
            season_code.append(''.join(month_codes[start[i] - 1 : 12] + month_codes[0 : end[i]]))
        
        
    return season_code
    
def calculate_translation(data):
    
    # Calculating translation
    data['translation'] = 6 - data['month']
    
    
    return data

def adjustMonth(month, translation):
    
    if month > 12:
        month = month - 12
    elif month < 1:
        month = month + 12
        
    return month
    
def apply_translation(data, trans_data):
    
    # Grouping the data by cell_ID and applying corresponding translation of month 
    # column
    new_data = []
    for cell, data_group in data.groupby(by = 'cell_ID'):
        # print(cell)
        # print(data_group)
        
        trans_group = trans_data.loc[trans_data['cell_ID'] == cell, 'translation'].item()
        
        # print(trans_group)
        
        data_group['month_trans'] = data_group['month'] + trans_group
        
        # print(data_group)
        # Adjusting month_trans label to be within 1 and 12
        data_group['month_trans'] = data_group['month_trans'].apply(adjustMonth, translation = trans_group)
        
        # print(data_group)
        
        # Appending modified group to list of dataframes
        new_data.append(data_group)
        
    
    # Converting list of lists to unqieu dataframe
    new_data = pd.concat(new_data, ignore_index = True)
    data = new_data
    
    # Freeing memory
    del new_data
    # print(data[data['cell_ID'] == 132])
    
    return data
    

    
    
def chmap_seasonLength(grid_size, temp, criteria):
    
    data_fp = '../data/processed/MBfogo_polygon_grid{}km_cerrado_2019_criteria0{}.csv'.format(grid_size, criteria)
    
    grid_fp = '../data/grids/cerrado_grid_{}km_cropped.shp'.format(grid_size)

    # Reading files
    data = pd.read_csv(data_fp)
    grid = gpd.read_file(grid_fp)
    
    # Set cell_ID column name depending on critera
    # First, get name based on criteria
    if criteria == 1:
        variable = 'cell_inters'
    elif criteria == 2:
        variable = 'cell_assign'
        
    # Then, change to 'cell_ID'
    data = data.rename(columns = {variable : 'cell_ID'})
    
    # Defining set of quantiles
    quantiles = [(0.25, 0.75, 50), (0.2, 0.8, 60), (0.15, 0.85, 70), (0.1, 0.9, 80)]
    # quantiles_labels = ['50', '60', '70', '80']
    
    # First getting information of month with more fires for each cell
    month_info = select_monthMoreFires_perCell(data)
    month_info = month_info.drop(columns = 'npol')
    month_info = calculate_translation(month_info)
    # print(month_info)
    
    data = apply_translation(data, month_info)
    
    # Calculate fire season parameters for each different set of quantiles
    for q in quantiles:
        
        # Start of fire season in translated label
        start_season = data[['cell_ID', 'month_trans']].groupby(by = 'cell_ID').quantile(q[0], interpolation = 'nearest').reset_index()
        
        # End of fire season in translated label
        end_season = data[['cell_ID', 'month_trans']].groupby(by = 'cell_ID').quantile(q[1], interpolation = 'nearest').reset_index()
        
        # print(start_season)
        # print(end_season)
        
        # Merging datasets
        season = pd.merge(start_season, end_season, how = 'outer', on = 'cell_ID',
                          suffixes = ('_start', '_end'))
        
        # Calculating season length
        season['length'] = season['month_trans_end'] - season['month_trans_start'] + 1
        
        # Merging to month_info
        season = pd.merge(season, month_info[['cell_ID', 'translation']], on = 'cell_ID')
        
        # Getting the start and end of the fires season in 'original month labels'
        season['month_start'] = season['month_trans_start'] - season['translation']
        season['month_end'] = season['month_trans_end'] - season['translation']
        
        season['month_start'] = [adjustMonth(season.loc[i, 'month_start'], season.loc[i, 'translation']) for i in range(len(season))]
        season['month_end'] = [adjustMonth(season.loc[i, 'month_end'], season.loc[i, 'translation']) for i in range(len(season))]
        
        
        # print(season)
        
        # Creating season code
        # print(season.tail(20))
        # print(season['month_end'].tail(20))
        
        
        season['code'] = season_code(season['month_start'], season['month_end'])
        
        # print(season)
        # print(season['month_trans_start'].unique())
        
        # Merge to grid
        grid_toPlot = grid.merge(season, left_on = 'id', right_on = 'cell_ID', how = 'outer')

        
        # Generating map for start of fire season
        chmap.chmap(grid_toPlot, 'month_start', 'Month', 
                title = 'Month of fire season start per cell in Cerrado in 2019.\n Using interquentile of {}%\n Grid of {} km and criteria {}'.format(q[2], grid_size, criteria),
                filename = 'MBfogo_seasonStart_map_iq{}_grid{}km_criteria{}.png'.format(q[2], grid_size, criteria),
                categorical = True, colmap = 'viridis')
        
        # Generating map for end of fire season
        chmap.chmap(grid_toPlot, 'month_end', 'Month', 
                title = 'Month of fire season end per cell in Cerrado in 2019.\n Using interquentile of {}%\n Grid of {} km and criteria {}'.format(q[2], grid_size, criteria),
                filename = 'MBfogo_seasonEnd_map_iq{}_grid{}km_criteria{}.png'.format(q[2], grid_size, criteria),
                categorical = True, colmap = 'viridis')
        
        # Generating map for length of fire season
        chmap.chmap(grid_toPlot, 'length', 'Number of months', 
                title = 'Length of fire season per cell in Cerrado in 2019.\n Using interquentile of {}%\n Grid of {} km and criteria {}'.format(q[2], grid_size, criteria),
                filename = 'MBfogo_seasonLength_map_iq{}_grid{}km_criteria{}.png'.format(q[2], grid_size, criteria),
                categorical = True, colmap = 'viridis')
        
        # Histogram for start of fire season
        hist.hist(grid_toPlot, 'month_start', 
            title = 'Histogram of number of cells by start of the fire season (Cerrado, 2019) \n Using interquentile of {}%\n Grid of {}km'.format(q[2], grid_size),
            bins = np.arange(1, 13+1) -0.5,
            #rangeb = (1, 11),
            xlabel = 'Start of the fire season (month)', ylabel = 'Number of cells',
            filename = 'MBfogo_seasonStart_hist_iq{}_grid{}km.png'.format(q[2], grid_size))        
        
        # Histogram for end of the fire season
        hist.hist(grid_toPlot, 'month_end', 
            title = 'Histogram of number of cells by end of the fire season (Cerrado, 2019) \n Using interquentile of {}%\n Grid of {}km'.format(q[2], grid_size),
            bins = np.arange(1, 13+1) -0.5,
            # rangeb = (1, 12),
            xlabel = 'End of the fire season (month)', ylabel = 'Number of cells',
            filename = 'MBfogo_seasonEnd_hist_iq{}_grid{}km.png'.format(q[2], grid_size))   
        
        # Histogram for length of fire season
        hist.hist(grid_toPlot, 'length', 
            title = 'Histogram of number of cells by length of the fire season (Cerrado, 2019) \n Using interquentile of {}%\n Grid of {}km'.format(q[2], grid_size),
            bins = np.arange(1, 13+1) -0.5,
            # rangeb = (1, 12),
            xlabel = 'Length of the fire season (number of months)', ylabel = 'Number of cells',
            filename = 'MBfogo_seasonLength_hist_iq{}_grid{}km.png'.format(q[2], grid_size))   
        
        # # Generating map for season type
        # chmap.chmap(grid_toPlot, 'code', 'Season code', 
        #         title = 'Fire season per cell in Cerrado in 2019.\n Using interquentile of {}%\n Grid of {} km and criteria {}'.format(q[2], grid_size, criteria),
        #         filename = 'MBfogo_seasonCode_map_iq{}_grid{}km_criteria{}.png'.format(q[2], grid_size, criteria),
        #         categorical = True, colmap = 'viridis')
        
        
def averageFreq(data, list_cols):
    
    # Initialising columns at 0
    data['mean_freq'] = 0
    
    # calculate weighted sum
    for col in list_cols:
        data['mean_freq'] += data[col] * int(col)
        
    # Finally, divide by number of pixels per cell to obtain the mean frequency
    data['mean_freq'] = data['mean_freq'] / data['npixels']
    
    return data

def quantileFreq(data, list_cols):
    
    # Create column with 0s
    data['95th_freq'] = 0
    

    # For each row, calculate the frequencies' 95th quantile
    for i in range(len(data.index)):
        
        # List of frequencies per pixel - intialise
        freq = []
        
        # For each column, add the corresponding frequency as many times as pixels with
        # such value in cell
        for col in list_cols:
            freq.extend([int(col)] * data.loc[i, col])
            
        # Now we have a list of frequencies, where each frequency appears as many
        # times as pixels with this frequency in the cell
        
        # Calculate 95th quantile of such a list
        data.loc[i, '95th_freq'] = np.quantile(freq, 0.95)
        
        
    return data
        
        
def chmap_freq(grid_size, temp):
    
    # File with the frequency data
    data_fp = '../data/processed/MBfogo_freq_grid{}km_cerrado_2019.csv'.format(grid_size)
        
    grid_fp = '../data/grids/cerrado_grid_{}km_clipped.shp'.format(grid_size)

    # Reading files
    data = pd.read_csv(data_fp)
    grid = gpd.read_file(grid_fp)  
    
    # First, create list of relevant column names (this is kind of manual after
    # checking max frequency in file)
    list_cols = [str(i) for i in range(36 + 1)]
    
    # Then, we calculate the total number of pixels
    data['npixels'] = data[list_cols].sum(axis = 1)
    
    # Calculate the average frequency
    data = averageFreq(data, list_cols)
    
    # Calculate the 95th quantile
    data = quantileFreq(data, list_cols)
    
    # Calculate area percentage per frequency value        
    # Calculate percentages
    data.loc[:, list_cols] = 100 * data.loc[:, list_cols].div(data['npixels'], axis = 0)
    
    # print(data[list_cols])
    
    # Creating columns for the maps that we want to plot: the sum of % for all frequency > x.
    # With x = 5, 10, 15, 20, 25, 30
    x = [5, 10, 15, 20, 25, 30]
    for i in x:
        data['freq_gt_{}'.format(i)] = data[list_cols[i: len(list_cols)]].sum(axis = 1)
        
        
    # Finally, merging relevant columns to grid
    grid = grid.merge(data[['id', '0', 'freq_gt_5', 'freq_gt_10', 'freq_gt_15', 'freq_gt_20', 'freq_gt_25', 'freq_gt_30', 'mean_freq', '95th_freq']], 
                      on = 'id', how = 'outer')
    
        
    # Now we have all the information to produce the maps.
    
    # 1. Percentage of area never burned.
    chmap.chmap(grid, '0', 'Percentage area', 
                title = 'Percentage of area never burned in Cerrado in 2019.\n Grid of {} km'.format(grid_size),
                filename = 'MBfogo_freq_neverBurned_map_grid{}km.png'.format(grid_size),
                categorical = False, colmap = 'viridis')
    
    # 2. Percentage of area that burned x times or more 
    for i in x:
        chmap.chmap(grid, 'freq_gt_{}'.format(i), 'Percentage area', 
                    title = 'Percentage of area burned {} times or more in Cerrado in 2019.\n Grid of {} km'.format(i, grid_size),
                    filename = 'MBfogo_freq_gt{}_map_grid{}km.png'.format(i, grid_size),
                    categorical = False, colmap = 'viridis')
        
    # 3. Average frequency
    # Choropleth map
    chmap.chmap(grid, 'mean_freq', 'Frequency', 
                title = 'Average fire frequency per cell in Cerrado in 2019.\n Grid of {} km'.format(grid_size),
                filename = 'MBfogo_freq_average_map_grid{}km.png'.format(grid_size),
                categorical = False, colmap = 'viridis')
    
    # # What if we erase those cells were frequency is 0
    # grid = grid[grid['mean_freq'] > 0]
    print(grid['mean_freq'].max() + 1)
    # Histogram plot
    hist.hist(grid, 'mean_freq', 
            title = 'Histrogram of average fire frequency in 36 years period per grid cell (Cerrado, 2019) \n Grid of {}km'.format(grid_size),
            bins = int(grid['mean_freq'].max() + 1),
            #rangeb = (1, grid['mean_freq'].max()),
            xlabel = 'Frequency (years)', ylabel = 'Number of cells',
            filename = 'MBfogo_freq_average_hist_grid{}km.png'.format(grid_size))
    
    # 4. The 95th quantile
    chmap.chmap(grid, '95th_freq', 'Frequency', 
                title = '95th quantile of fire frequency per cell in Cerrado in 2019.\n Grid of {} km'.format(grid_size),
                filename = 'MBfogo_freq_q95_map_grid{}km.png'.format(grid_size),
                categorical = False, colmap = 'viridis')
    # TODO: set the number of bins to be the same as the max value in series minus minimum so one bins per category
    print(grid['95th_freq'].max() + 1)
    # Histogram plot
    hist.hist(grid, '95th_freq', 
            title = 'Histrogram of 95th quantile fire frequency in 36 years period per grid cell (Cerrado, 2019) \n Grid of {}km'.format(grid_size),
            bins = int(grid['95th_freq'].max() + 1),
            #rangeb = (1, grid['95th_freq'].max()),
            xlabel = '95th quantile of fire frequency (years)', ylabel = 'Number of cells',
            filename = 'MBfogo_freq_q95s_hist_grid{}km.png'.format(grid_size))
    
    
def chmap_landuse(grid_size):
    
    # File with the frequency data
    data_fp = '../data/processed/MBfogo_landuse_grid{}km_cerrado_2019.csv'.format(grid_size)
        
    grid_fp = '../data/grids/cerrado_grid_{}km_cropped.shp'.format(grid_size)

    # Reading files
    data = pd.read_csv(data_fp)
    grid = gpd.read_file(grid_fp)
    
    # First, calculate percentages
    data.iloc[:, data.columns.get_loc('1') : data.columns.get_loc('49')] = \
        100 * data.iloc[:, data.columns.get_loc('1') : data.columns.get_loc('49')].div(data['npixels'], axis = 0)
        
    # Add up total percentage for anthropic categories
    data['anthropic'] = data[['14', '24', '30']].sum(axis = 1)
    
    
    # Columns for which we want to produce a map
    categories = ['anthropic', '4', '3', '12', '11', '15', '18']
    
    # Getting dictionaries of MapBiomas la =nd use legend codes
    legend_codes = mblegend.legend_names_MBc6()
    
    # Labels for each category/group of categories
    labels = ['Anthropic'] + [legend_codes[i] for i in ['4', '3', '12', '11', '15', '18']]
    
    # Finally, merging relevant columns to grid
    grid = grid.merge(data[['id'] + categories], 
                      on = 'id', how = 'outer')
    
    # Produce map per category
    for i in range(len(categories)):
        chmap.chmap(grid, categories[i], 'Percentage area', 
                    title = 'Area percentage of {} land-use in Cerrado in 2019.\n Grid of {} km'.format(labels[i], grid_size),
                    filename = 'MBfogo_landuse_{}_map_grid{}km.png'.format(labels[i], grid_size),
                    categorical = False, colmap = 'viridis')
        
    # Histogram plot of anthropic cover
    hist.hist(grid, categories[0], 
            title = 'Histrogram of percentage anthropic area per grid cell (Cerrado, 2019) \n Grid of {}km'.format(grid_size),
            bins = 20,
            xlabel = 'Percentage area (%)', ylabel = 'Number of cells',
            filename = 'MBfogo_landuse_anthropic_hist_grid{}km.png'.format(grid_size))
        

def chmap_grid(variable, grid_size, temp, criteria = None):
    
    # Call appropriate function depending on variable
    if variable == 'burnedArea':
        chmap_burnedArea(grid_size, temp)
        
    elif variable == 'numPolygons':
        chmap_numPolygons(grid_size, temp, criteria)
        
    elif variable == 'sizeQuantile':
        chmap_sizeQuantile(grid_size, temp, criteria)
        
    elif variable == 'seasonality_moreFires':
        chmap_seasonality_moreFires(grid_size, temp, criteria)
    
    elif variable == 'seasonality_largestFire':
        chmap_seasonality_largestFire(grid_size, temp, criteria)
    
    elif variable == 'season_length':
        chmap_seasonLength(grid_size, temp, criteria)
        
    elif variable == 'frequency':
        chmap_freq(grid_size, temp)
    
    elif variable == 'landuse':
        chmap_landuse(grid_size)
    
    


"""
MAIN
"""
def main():
    chmap_grid('sizeQuantile', 50, 'y', criteria = 2)
    
    # data_fp = '../data/processed/MBfogo_polygon_grid50km_cerrado_2019_criteria02.csv'
    # data = pd.read_csv(data_fp)
    # print(data[data['area'] == data['area'].max()])
    
if __name__ == '__main__':
    main()

