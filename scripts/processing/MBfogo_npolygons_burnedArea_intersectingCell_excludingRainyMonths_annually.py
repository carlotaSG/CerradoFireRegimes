# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 17:10:48 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that counts the number of polygons assigned to a cell not falling in the 
peak of the rainy season (the three months with highest average rain in the 
period). As well, it get the burned area per month (from pixel method) file, 
casts the rainy season months as 0, combines this information with the above and
outputs to file.

Input files:
    - ../../data/processed/summary_tables/<nyears>_year_periods/chirps_meanMonthlyPrec_period<>.csv
    - ../../data/processed/MBfogo_c10_cerrado_polygons_b150_CM/MBfogo_c10_cerrado_polygons_<>.csv
    - ../../data/processed/summary_tables/<nyears>_year_periods/burnedArea_fireCount_intersection_monthFromPixel_period<>.csv
    
Output files:
    - ../../data/processed/summary_tables/<nyears>_years_period/burnedArea_monthFromPixel_fireCount_assignedCell_exclude<>RainyMonths_grid<>_period<>.csv
"""

"""
Imports
"""
import pandas as pd


"""
Functions
"""


"""
MAIN
"""
def main():
    
    # ----------------- User inputs -------------------------------------------
    
    # UPDATE (25/10/2023):
    # Allowing user to choose MapBiomas FOGO collection
    collection = 2
    
    in_dir = '../../data/processed/'
    
    prec_data_fp = in_dir + 'climate_data/chirps_grid{}_{}.csv'
    
    polygons_fp = in_dir + 'MBfogo_c{}0_cerrado_polygons_b150_CM/MBfogo_c{}0_cerrado_polygons_{}.csv'.format(collection, collection, '{}')
    
    ba_data_fp = in_dir + 'summary_tables/annual/burnedArea_fireCount_intersection_monthFromPixel_grid{}_year{}.csv'
    
    output_fp = in_dir + 'summary_tables/annual/burnedArea_monthFromPixel_fireCount_assignedCell_exclude{}RainyMonths_grid{}_year{}.csv'
    
    

    # # Number of periods
    # n_periods = 4
    
    # if n_periods == 2:
    #     # Year periods to work on 
    #     periods = [(1985, 2002), (2003, 2020)]
        
    #     # Number of years per period
    #     n_years = 18
        
    # elif n_periods == 4:
    #     # Year periods to work on 
    #     periods = [(1985,1993), (1994, 2002), (2003, 2011), (2012, 2020)]
        
    #     # Number of years per period
    #     n_years = 9
    
    nmonths = 3
    
    min_polSize = 0.03
    
    grid_id = '30km'
    
    start_year, end_year = 1985, 2020
    
    grid_colname = 'cellID_{}'.format(grid_id)
    
    # Working one year at a time
    for year in range(start_year, end_year + 1):
        
        print('Working on year {}'.format(year))
        
        # Read average precipitation per month and cell for this year
        prec_data = pd.read_csv(prec_data_fp.format(grid_id, year))
        
        # This file already only contains data for those cells with more than 
        # a certain percentage of their area inside the Cerrado
        
        # Getting the list of cells
        list_cells = prec_data['id'].unique()
        
        cells_rainy_months = {}
        for cell in list_cells:
            
            # Find the x months with more average precipitation
            rainy_months = prec_data.loc[prec_data['id'] == cell, ['month', 'mean_prec']]
            rainy_months = rainy_months.sort_values(by = 'mean_prec', ascending = False).reset_index(drop = True)
            rainy_months = rainy_months.loc[:nmonths-1, 'month'].tolist()
            
            cells_rainy_months[cell] = rainy_months
            
        
        # Reading polygons for this year
        # polygons = []
        # for year in range(periods[i-1][0], periods[i-1][1] + 1):
            
        polygons_year = pd.read_csv(polygons_fp.format(year))
            
        # Filtering out polygons smaller than threshold
        polygons_year = polygons_year[polygons_year['area'] >= min_polSize]
            
        # UPDATE (25/10/2023):
        # With the new grids, we do not exclude the border polygons
        # # Filtering out border polygons
        # polygons_year = polygons_year[polygons_year['border_pol'] == 0]
        
        # Selecting columns of interest
        polygons_year = polygons_year[['month', grid_colname]]
            
        # Renaming grid columns
        polygons_year = polygons_year.rename(columns = {grid_colname: 'id'})
            
        # Filter out cells we are not interested in
        polygons_year = polygons_year[polygons_year['id'].isin(list_cells)]
            
            # polygons.append(polygons_year)
            
            # del polygons_year

        # polygons = pd.concat(polygons, ignore_index = True)      
        
        # Counting polygons per cell and month in the year
        polygons_year = polygons_year.groupby(by = ['id','month']).size().reset_index(name = 'npolygons')
        polygons_year = polygons_year.sort_values(['id','month']).reset_index(drop = True)
        
        # Adding column with year and reordering
        polygons_year['year'] = year
        polygons_year = polygons_year[['id', 'year', 'month', 'npolygons']]
        
        
        # Read burned area per cell, year and month
        ba_data = pd.read_csv(ba_data_fp.format(grid_id, year))
        
        # Filter out cells we are not interested in
        ba_data = ba_data[ba_data['id'].isin(list_cells)]
        
        # Dropping unnecessary columns
        ba_data = ba_data.drop(columns = ['npolygons','area_T'])
        
        # Casting from wide to long format
        ba_data = pd.melt(ba_data, 
                          id_vars = ['id','year'], 
                          var_name = 'month', value_name = 'area')
        ba_data['month'] = ba_data['month'].apply(int)
        
        # Sorting values
        ba_data = ba_data.sort_values(['id','year','month']).reset_index(drop = True)
        
        
        # Merge burned area and polygon data
        ba_data = ba_data.merge(polygons_year, on = ['id','year','month'], how = 'left')
        
        # Filling Nas with 0s
        ba_data = ba_data.fillna(0)
        
        # Casting data for the wettest months as 0 for each cell
        for cell in list_cells:
            
            ba_data.loc[(ba_data['id'] == cell) & (ba_data['month'].isin(cells_rainy_months[cell])), ['area','npolygons']] = 0
        
        
        print('...saving file...\n')
        ba_data.to_csv(output_fp.format(nmonths, grid_id, year), index = False)
        

    
    
    
    
    
    
if __name__ == '__main__':
    main()

