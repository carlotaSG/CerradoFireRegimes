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
    grid_id = '30km'
    
    # Number of periods we are working with
    n_periods = 4#2
    
    # Number of years per period
    n_years = 9#18
    
    # Threshold area
    area_thresh = 0.0
    
    # # Extra cells to include and cells to exclude
    # include_cells = [567, 568]
    # exclude_cells = [1451]
    
    # Data files and the columns we want from every file
    # data_fp = '../../data_new_grids/processed/summary_tables/annual/population_density_grid{}_{}.csv'.format(grid_id, '{}')
    # data_cols = ['areaw_popdensity'] 
    
    # data_fp = '../../data_new_grids/processed/summary_tables/annual/livestock_density_grid{}_{}.csv'.format(grid_id, '{}')
    # data_cols = ['areaw_dens']
    
    # data_fp = '../../data_new_grids/processed/landscape_metrics/MBlulc_landscapeMetrics_grid{}_{}.csv'.format(grid_id, '{}')
    # data_cols = ['total_area', 'number_of_patches', 'patch_density', 'total_edge', 'edge_density', 'landscape_shape_index', 'area_mn', 'area_cv']
    
    data_fp = '../../data_new_grids/processed/landscape_metrics/MBlulc_landscapeCellMetrics_grid{}_{}.csv'.format(grid_id, '{}')
    data_cols = ['class_val', 'total_area', 'proportion_of_landscape', 'number_of_patches', 'patch_density', 
                  'largest_patch_index', 'total_edge', 'edge_density', 
                  'landscape_shape_index', 'area_mn', 'area_md', 'area_cv']
    
    # data_fp = '../../data_new_grids/processed/protected_areas/protected_areas_{}_{}.csv'.format(grid_id, '{}')
    # data_cols = ['pa_total']
    
    # output_fp = '../../data_new_grids/processed/summary_tables/{}_year_periods/population_density_allPeriods.csv'.format(n_years)
    # output_fp = '../../data_new_grids/processed/summary_tables/{}_year_periods/livestock_density_allPeriods.csv'.format(n_years)
    # output_fp = '../../data_new_grids/processed/summary_tables/{}_year_periods/landscape_metrics_grid{}_allPeriods.csv'.format(n_years, grid_id)
    output_fp = '../../data_new_grids/processed/summary_tables/{}_year_periods/landscape_cell_metrics_grid{}_allPeriods.csv'.format(n_years, grid_id)
    # output_fp = '../../data_new_grids/processed/summary_tables/{}_year_periods/protected_areas_allPeriods.csv'.format(n_years)
    
    landsCellMetrics = True
    
    
    # -------------------- Reading and formatting data -----------------------
    
    # Reading and formatting the data
    data = readFormat_dataYears(data_fp, start_year, end_year, grid_id, data_cols,
                             periods = True, nperiods = n_periods, month = False,
                             drop_geometry = True,
                             thresh = area_thresh)
    data = data.drop(columns = '%area_crop').reset_index(drop = True)
    # print(data)
    
    # ----------------- Summarising the data ---------------------------------
    
    # Variable to store the cell's data
    data_out = []
    
    if landsCellMetrics == True:
        cellMetrics = ['class_val']
        col_position = 3
    else:
        cellMetrics = []
        col_position = 2
    
    # Calculate the average and the change in variable
    for (cell, group) in tqdm(data.drop(columns = 'year').groupby(by = ['id','period'] + cellMetrics)):
        # print(cell)
        # print(group.reset_index(drop = True))
    
        group = group.reset_index(drop = True)
        
        # Calculating the average data values in period
        part1 = group.drop(columns = ['id', 'period'] + cellMetrics).apply('median').tolist()
        
        # Calculate the change in data value in period
        part2 = group.drop(columns = ['id', 'period'] + cellMetrics).apply(lambda x: x[x.tail(1).index].item() - x[0]).tolist()
        
        # print(part1)
        
        # print(part2)
        
        # Setting the new column names as the ones we had (for the averages) 
        # and the suffix _delta for the changes in value over the period
        new_cols = list( group.columns) + [el + '_delta' for el in group.columns[col_position:]]
        # print(new_cols)
        # print(cell)
        # print([cell[0], cell[1]] + part1 + part2)
        
        # Casting as a dataframe
        if landsCellMetrics == False:
            cell_out = pd.DataFrame(pd.Series([cell[0], cell[1]] + part1 + part2, index = new_cols)).transpose()
        elif landsCellMetrics == True:
            cell_out = pd.DataFrame(pd.Series([cell[0], cell[1], cell[2]] + part1 + part2, index = new_cols)).transpose()
        # Accumulating in dataframe
        data_out.append(cell_out)
        # break 
    
    # Casting list of dataframes as unique dataframe
    data_out = pd.concat(data_out, ignore_index = True)
    
    # Saving file
    data_out.to_csv(output_fp, index = False)

    # print(data.drop(columns = 'year').groupby(by = ['id','period']).apply(lambda x: x[x.tail(1).index].item() - x[0].item()))
    
    
if __name__ == '__main__':
    main()