# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 10:29:46 2023

@author: scat8298
"""

"""
Imports
"""
import pandas as pd
import sys
sys.path.append('../utils/')
from readFormat_dataYears import readFormat_dataYears
from tqdm import tqdm


"""
Functions
"""

"""
MAIN
"""
def main():
    
    # ------------------------- User inputs ---------------------------------
    
    # Year period
    start_year = 1985
    end_year = 2020
    
    # Grid to use
    grid_id = '50km'
    
    # Number of periods we are working with
    n_periods = 4
    
    # Number of years per period
    n_years = 9
    
    # Threshold area
    area_thresh = 0.0
    
    # # Extra cells to include and cells to exclude
    # include_cells = [567, 568]
    # exclude_cells = [1451]
    
    # Data files and the columns we want from every file
    lulc_fp = '../../data_new_grids/processed/lulc_data/MBlanduse_percArea_subsetLULC_grid{}_{}.csv'.format(grid_id, '{}')
    lulc_cols = ['natural', 'anthropic', '4', '3', '12', '11', '15', '18', '9', '21']
    
    output_fp = '../../data_new_grids/processed/summary_tables/{}_year_periods/lulc_percArea_grid{}_allPeriods.csv'.format(n_years, grid_id)
    
    # -------------------- Reading and formatting data -----------------------
    
    # Reading and formatting the data
    lulc = readFormat_dataYears(lulc_fp, start_year, end_year, grid_id, lulc_cols,
                             periods = True, nperiods = n_periods, month = False,
                             drop_geometry = True,
                             thresh = area_thresh)
    lulc = lulc.drop(columns = '%area_crop').reset_index(drop = True)
    # print(lulc)
    
    # ----------------- Summarising the data ---------------------------------
    
    # Variable to store the cell's data
    data_out = []
    
    # Calculate the average and the change in variable
    for (cell, group) in tqdm(lulc.drop(columns = 'year').groupby(by = ['id','period'])):
        # print(cell)
        # print(group.reset_index(drop = True))
    
        group = group.reset_index(drop = True)
        
        # Calculating the average LULC values in period
        part1 = group.drop(columns = ['id', 'period']).apply('mean').tolist()
        
        # Calculate the change in LULC value in period
        part2 = group.drop(columns = ['id', 'period']).apply(lambda x: x[x.tail(1).index].item() - x[0]).tolist()
        
        # print(part1)
        
        # print(type(group.columns))
        
        # Setting hte new column names as the ones we had (for the averages) 
        # and the suffix _delta for the changes in value over the period
        new_cols = list( group.columns) + [el + '_delta' for el in group.columns[2:]]
        # print(new_cols)
        # print(cell)
        # print([cell[0], cell[1]] + part1 + part2)
        
        # Casting as a dataframe
        cell_out = pd.DataFrame(pd.Series([cell[0], cell[1]] + part1 + part2, index = new_cols)).transpose()
        
        # Accumulating in dataframe
        data_out.append(cell_out)
        # break 
    
    # Casting list of dataframes as unique dataframe
    data_out = pd.concat(data_out, ignore_index = True)
    
    # Saving file
    data_out.to_csv(output_fp, index = False)

    # print(lulc.drop(columns = 'year').groupby(by = ['id','period']).apply(lambda x: x[x.tail(1).index].item() - x[0].item()))
    
    
if __name__ == '__main__':
    main()