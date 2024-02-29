# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 12:46:18 2023

@author: scat8298

Script that counts the number of fire polygons assigned to each cell in the 
period, as well as adding up the total area of the polygons assigned to the cell
in the period.

Input files:
    - ../../data/processed/MBfogo_c10_cerrado_polygons_b150_CM/MBfogo_c10_cerrado_polygons_{}.csv

Output files:
    - ../../data/processed/summary_tables/{}_year_periods/burnedArea_fireCount_assigned_grid{}_period{}.csv
"""

"""
Imports
"""
import pandas as pd

"""
MAIN
"""
    
def main():
    
    
    #----------------- User inputs -------------------------------------------
    
    # UPDATE 25/10/2023
    # Choosing MapBiomas FOGO collection
    collection = 2
    
    # Filepaths to polygons and their associated data 
    polygons_data_fp = '../../data/processed/MBfogo_c{}0_cerrado_polygons_b150_CM/MBfogo_c{}0_cerrado_polygons_{}.csv'.format(collection, collection, '{}')
    
    # Number of periods
    n_periods = 4
    
    if n_periods == 2:
        # Groups of years - script will generate a summary file per group
        group_years = [(1985, 2002), (2003, 2020)]
        
        # Number of years per period
        n_years = 18
        
    elif n_periods == 4:
        # Groups of years - script will generate a summary file per group 
        group_years = [(1985,1993), (1994, 2002), (2003, 2011), (2012, 2020)]
        
        # Number of years per period
        n_years = 9
    
    grid_id = '30km'
    
    # Output filepath
    output_fp = '../../data/processed/summary_tables/{}_year_periods/burnedArea_fireCount_assigned_grid{}_period{}.csv'
    
    annual_output_fp = '../../data/processed/summary_tables/annual/burnedArea_fireCount_assigned_grid{}_year{}.csv'
    
    #-------------- Area calculation -----------------------------------------

    # Counter to loop over the periods
    counter = 1
    
    # Looping over groups of years
    # We work on a group of years at a time
    for (start_year, end_year) in group_years:
        
        period = '({}, {})'.format(start_year, end_year)
        
        print('Working on group ' + period)
        
        # Creating empty list to accumulate the burned area and count of polygons 
        # for each cell in each year and month
        df_period = []
        
    
        # For each year in group, calculate the burned area inside each cell
        for year in range(start_year, end_year + 1):
            print('')
            print('Working on year {}'.format(year))

            # The data file for this year
            df_year = pd.read_csv(polygons_data_fp.format(year))
            
            # Filter to have only polygons above a certain are threshold
            df_year = df_year[df_year['area'] >= 0.03]
            
            # Select only columns of interest
            df_year = df_year[['cellID_{}'.format(grid_id), 'month', 'area']]
            
            # Rename cell id columns
            df_year = df_year.rename(columns = {'cellID_{}'.format(grid_id): 'id'})
            
            # Counting number of polygons per cell
            npolygons_cell = df_year[['id']].groupby(by = 'id').size().reset_index(name = 'npolygons')
            
            # Adding up the area for all cell instances for the same month
            df_year = df_year.groupby(by = ['id', 'month']).sum().reset_index()
            
            # print(df_year)
            
            # Convert month to string
            df_year['month'] = df_year['month'].apply(str)
            
            # Convert month column to columns (long-to-wide) with area as the content
            df_year = df_year.pivot(index = 'id', columns = 'month', values = 'area').reset_index()
            
            # print(df_year)
            
            # Filling nas with 0
            df_year = df_year.fillna(0)
            
            # Check if missing month columns, therwise add
            complete_cols = set(['id'] + [str(m) for m in range(1, 12+1)])
            actual_cols = set(df_year.columns)
            # Missing columns
            missing_cols = list(sorted(complete_cols - actual_cols))
            
            if len(missing_cols) > 0:
                df_year[missing_cols] = 0

            
            # Calculate total area per cell
            df_year['area_T'] = df_year[[str(m) for m in range(1, 12 + 1)]].sum(axis = 1)
            
            # Merge with the number of polygons per cell
            df_year = df_year.merge(npolygons_cell, on = 'id', how = 'outer')
            
            # -----------------------------------------------------------------
            # UPDATE 17/11/2023: I need the annual information per cell!!
            
            # Create column with year id
            df_year['year'] = year
            # Casting columns in order
            df_year = df_year[['id', 'year', 'npolygons','area_T'] + [str(m) for m in range(1, 12+1)]]
            # Saving to file
            df_year.to_csv(annual_output_fp.format(grid_id, year), index = False)
            # -----------------------------------------------------------------
            
            # making sure all columns are in the same order
            df_year = df_year[['id', 'npolygons','area_T'] + [str(m) for m in range(1, 12+1)]]
            
            # Accumulating this information for the period
            df_period.append(df_year)
            
            # Freeing memory
            del df_year, npolygons_cell
            
        # Converting list of dataframes to dataframes
        df_period = pd.concat(df_period, ignore_index = True)
    
        # Total area and number of polygons per cell and category
        df_period = df_period.groupby(by = 'id').sum().reset_index()
        
        # Now we add a column to each indicating the period
        df_period['period'] = period
        
        # Adding column to indicate the grid 
        df_period['grid_id'] = grid_id

        
        # Reshuffling columns
        df_period = df_period[['id','period','npolygons','area_T'] + [str(m) for m in range(1, 12+1)]]
        
        print('Saving data for all years in period...\n')
        
        # Saving dataframes
        df_period.to_csv(output_fp.format(n_years, grid_id, counter), index = False)
        
        counter += 1
            
            
        

if __name__ == '__main__':
    main()

