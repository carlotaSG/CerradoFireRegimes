# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 12:17:37 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculates the fire return interval as the inverse of the average 
annual percentage of the flammable area burned in a period. Then, the frequency
is calculated as the inverse of this quantity.

Input files:
    - ../../data/processed/lulc_data/MBlanduse_grid50km_cerrado_{}.csv
    - ../../data/processed/MBfogo_c10_cerrado_polygons_b150_CM/MBfogo_c10_cerrado_polygons_{}.{}
    - ../../data/raw/MBfogo_ba/MBfogo_c10_cerrado_{}.tif

Output files:
    - ../../data/processed/summary_tables/{}_years_period/fireFrequency_returnInterval_period{}.csv    


"""

"""
Imports
"""
import geopandas as gpd
import pandas as pd
from multiprocessing import Pool
from pytictoc import TicToc
t = TicToc()

import sys
sys.path.append('../utils')

import assignRasterValue_toPolygon
import parallel_prep as pprep


"""
Functions
"""           
        
def countCells_area_intersectingPolygon(slices, polygon_fp, polygon_data_fp, col_names, raster_fp, grid):
    """
    Script that calculates the sum of the areas of polygons intersecting each 
    cell. It works on a set  (slice) of polygons at a time: it reads the slice of
    polygons, finds their intersection with each cell, calculates the area of 
    the intersection, and returns a datafrmae where each row is a polygon 
    intersection, identified with the cell it is found in.

    Parameters
    ----------
    slices : List of integers
        Integers indicating the rows of hte polygons to read from.
    polygon_fp : string
        Filepath to the shapefile containing the polygon geometries.
    polygon_data_fp : string
        Filepath to the CSV containing the polygon data.
    col_names : list of strings
        List of columns to read from the polygon_dat_fp file.
    raster_fp : string
        Filepath to the raster file of burned area.
    grid : GeoDataFrame
        The grid of cells against which we'll check what cells intersect with 
        each polygon.

    Returns
    -------
    df_area_cell_monthPixel : DataFrame
        Dataframe containing the area of each polygon intersection with each cell.

    """
    
    
    # Read chunk of polygons in shapefile
    polygons = gpd.read_file(polygon_fp, rows = slice(slices[0], slices[1]))
    
    # Read chunk of polygon data in csv file
    polygons_data = pd.read_csv(polygon_data_fp, names = col_names, skiprows = lambda x: x not in range(slices[0] + 1, slices[1]+1+1))
    # print(polygons_data)
    # Convert grid CRS to that of polygon
    grid = grid.to_crs(polygons.crs)
    
    # Creating empty lists to accumulate each polygon's area data per cell
    df_area_cell_monthPixel = []
    
    # Loop over polygons to assign them to cells individually
    for i in range(len(polygons)): 
        
        # Selecting only the polygon geometry, id, year and month
        pol_geom = polygons.loc[i, 'geometry']
        pol_month = polygons.loc[i,'month']
        pol_year = polygons.loc[i,'year']
        pol__id = polygons.loc[i, 'ID']

        # Select only polygons larger than 3 ha
        pol_area = polygons_data.loc[i, 'area']        
        if pol_area >= 0.03:
            
            # First, get the list of cells intersecting the polygon
            cell_inters = grid.loc[grid['geometry'].intersects(pol_geom), 'id'].apply(int).reset_index(drop = True)

            # UPDATE (25/10/2023)
            # Now we may have polygons that do not intersect any cells
            # In that case, we do not operate on them
            if len(cell_inters) > 0:
                # Calculate the geometries of each cell intersection
                cell_assign_geom = grid.loc[grid['id'].isin(cell_inters), 'geometry'].intersection(pol_geom).reset_index(drop = True)
                
                    
                # Then, we calculate the area of each intersection
                # For this we use the following function, which counts the number of
                # pixels inside each geometry per pixel category (i.e. month)
                # TODO: This returns a dataframe with columns: 1 | 2 | 3 | ... | 12 
                try:
                    exception_triggered = False
                    pixels_per_geom_month = assignRasterValue_toPolygon.countRasterValues_inListPolygons(
                        (cell_assign_geom, raster_fp, 1, 12, 0))
                except Exception as e:
                    with open("./exceptionFRI_30km.txt", "a") as myfile:
                        myfile.write(pol__id + '\n')
                        exception_triggered = True
                
                if exception_triggered == False:
    
                    # Converting number of pixels to area
                    area_per_geom_month = pixels_per_geom_month * 30 * 30 * 1e-06
                    
                    # Freeing memory
                    del cell_assign_geom, pixels_per_geom_month
                        
                    # Calculate total area (sum of the area per month)
                    # I just sum across all columns to have the area per geometry
                    area_per_geom = area_per_geom_month.sum(axis = 1)
                    
                    
        
                    # First, the dataframe where month area comes form pixel level information
                    # List of cells that intersect the polygon, along with the total intersection data (all months summed up)
                    area_cell_monthPixel = pd.DataFrame(
                        {'id': cell_inters,
                         'year': pol_year,
                         'area_T': area_per_geom})
                                 
                    # Appending to list of dataframes
                    df_area_cell_monthPixel.append(area_cell_monthPixel)
        
                    
                    # Freeing memory
                    del area_cell_monthPixel, 
                    area_per_geom, area_per_geom_month, cell_inters
            

                
    # If no polygons were larger than the threshold, , or none were intersecting cells
    if len(df_area_cell_monthPixel) > 0:
        # Cast list of dataframes as a unique dataframe for each one of the two options
        df_area_cell_monthPixel = pd.concat(df_area_cell_monthPixel, ignore_index = True)

    
    # Otherwise, we'll return empty dataframes
    if len(df_area_cell_monthPixel) == 0:
        df_area_cell_monthPixel = pd.DataFrame(columns = ['id', 'year', 'area_T'])
        
    
    return df_area_cell_monthPixel


"""
MAIN
"""

def main():
    
    # ----------------------- User inputs -------------------------------------
    
    # Grid size
    grid_size_units = '30km'
    
    
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
    
    # Number of processes to parallelise over
    num_process = 38
    
    lulc_area_fp = '../../data/processed/lulc_data/MBlanduse_grid{}_cerrado_{}.csv'.format(grid_size_units, '{}')
    
    output_fp = '../../data/processed/summary_tables/{}_year_periods/fireFrequency_returnInterval_grid{}_period{}.csv'
    
    # UPDATE (24/10/2023)
    # MapBiomas FOGO collection
    collection = 2
    
    
    # ------ Classifying polygon by cell for each year in period -------------- 
    
    
    # Grid filepath
    grid_fp = '../../data/shapes/cerrado_grid_{}_cropped.shp'.format(grid_size_units)
    # Reading grid
    grid = gpd.read_file(grid_fp)
    # Selecting only relevant columns
    grid = grid[['id', '%area_crop', 'withinCerr','geometry']]
    

    # Counter variable
    counter = 4# 1
    
    # We work on a group of years at a time
    for (start_year, end_year) in [group_years[3]]:
        t.tic()
        
        period = '({}, {})'.format(start_year, end_year)
        
        print('Working on group ' + period)
        
        # Creating empty lists where we accumulate the data for all years in group
        df_area_year_period = []
        
    
        # For each year in group, calculate the number of polygons intersecting 
        # each cell and the total burned area of the polygons (only the area
        # inside the cell)
        
        for year in range(start_year, end_year + 1):
            t.tic()
            print('')
            print('Working on year {}'.format(year))
        
            # ---------------------------------------------------------------------
            # Getting paths to input files
        
            # Filepath to files
            files = '../../data/processed/MBfogo_c{}0_cerrado_polygons_b150_CM/MBfogo_c{}0_cerrado_polygons_{}.{}'.format(collection, collection, '{}', '{}')
            # Filepath to fire polygon shapefile
            polygon_fp =  files.format(year, 'shp')
            
            # Filepath to fire polygon data
            polygon_data_fp = files.format(year, 'csv')
            
            # Filepath to MBfogo raster
            raster_fp = '../../data/raw/MBfogo_ba/MBfogo_c{}0_cerrado_{}.tif'.format(collection, year)
            
        
            # ---------------------------------------------------------------------
            # Assigning each polygon to a cell
        
            # We identify the cells that a polygon is intersecting with and the 
            # area of the polygon that falls within each cell by parallelising over groups 
            # of polygons. To achieve this, we first calculate the groups of polygons 
            # as slices of indices (row names). Each process in Pool will be in 
            # charge of a certain group of rows of the polygons file.
            
            print('Creating pool input...')
            # Reading the polygon's csv file to get the total number of rows.
            polygon_data = pd.read_csv(polygon_data_fp)

            # Get the column names
            polygon_data_colNames = polygon_data.columns
            

            # Only the total number of rows is needed to be able to produce partitions
            npolygons = int(len(polygon_data.index))

            del polygon_data
            
            # Generate the slices
            slices = pprep.list_slices(npolygons, num_process*4)
            

            # Generating input to Pool using the slices
            pool_input = list(zip(slices, [polygon_fp]*len(slices), 
                                  [polygon_data_fp]*len(slices), [polygon_data_colNames]*len(slices), 
                                  [raster_fp]*len(slices), 
                                  [grid]*len(slices)))
    

            with Pool(num_process) as p: 
                
                print("Starting parallel computation in countCells_area_intersectingPolygon...")
                
                # Assign polygon to a certain cell, return two dataframes
                # Left: the area per month is calculated from the pixels in polygon
                # Right: the area per month is just the total area of the polygon in cell
                #        assigned to the polygon's month
                df_area_year = p.starmap(countCells_area_intersectingPolygon, pool_input)


                print('Finished parallel computation.')
                
                
            # Converting list of dataframes to a single dataframe
            df_area_year = pd.concat(df_area_year, ignore_index = True)
            
            # These dataframes contain the total area per cell and year
            # But cells may be duplicated in the dataframe (polygons in different slices intersected the same cell)
            # So, we get sum to get the total count and area
            df_area_year = df_area_year.groupby(by = ['id', 'year']).sum().reset_index()
            
            
            # Reading the flammable data for the year
            flammable_area_year = pd.read_csv(lulc_area_fp.format(year))
            
            # Adding up the flammable data
            flammable_area_year['flammable'] = flammable_area_year[['1','11','12','13','14']].sum(axis = 1)
            # Subsetting relevant columns
            flammable_area_year = flammable_area_year[['id','flammable']]
            
            # Merge this data with the annual burned area
            df_area_year = df_area_year.merge(flammable_area_year, on = ['id'], how = 'left')
            
            # Calculate percentage of flammable area burned
            df_area_year['perc_ba'] = df_area_year['area_T'] / df_area_year['flammable']
            
            # Now we have all the data for each cell in year, we accumulate it to next year
            df_area_year_period.append(df_area_year)

            
        
            
            # Freeing memory
            del df_area_year
            
            t.toc('Year {} done in'.format(year))
            print('------------------------------------------------------------')
            
        print('Calculating fire return interval of all years in period...')
            
        # Now we have the data for all years in group as two lists of dataframes
        # First we convert them to dataframes
        df_area_year_period = pd.concat(df_area_year_period, ignore_index = True)
        
        df_area_year_period = df_area_year_period.sort_values(by = ['id','year']).reset_index(drop = True)


        # Calculate average percentage of flammable burned area (and total burned area, and average flammable area)
        df_area_year_period = df_area_year_period.groupby(by = 'id').agg({'area_T': 'mean',
                                                                         'flammable':'mean',
                                                                         'perc_ba':'mean'}).reset_index()
        # Calculate return interval
        df_area_year_period['fri'] = df_area_year_period['perc_ba'].apply(lambda x: 1/x)
        
        # Calculate frequency from fri
        df_area_year_period['fri_freq'] = df_area_year_period['fri'].apply(lambda x: n_years/x)
        
        # Adding column with period
        df_area_year_period['period'] = str(group_years[counter-1])
        
        # Reorder columns
        df_area_year_period = df_area_year_period[['id', 'period', 'area_T','flammable','perc_ba', 'fri', 'fri_freq']]
        
        # Save data to csv
        df_area_year_period.to_csv(output_fp.format(n_years, grid_size_units, counter), index = False)
        
        counter += 1
        # print(df_area_year_period.head(40))
        t.toc('Period calculated in')

        print('--------------------------------------------------------------')
        print('\n\n')
        
    

        
        
if __name__ == '__main__':
    main()

