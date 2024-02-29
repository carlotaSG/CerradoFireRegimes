# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 12:17:37 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculates the number of polygons intersecting each cell in each 
period, the total area of the intersections of all these polygons, as well as
the area of the intersections per month.

Functions:
    * countCells_area_intersectingPolygon - a version of this function is used 
                                            in script ./MBfogo_fireReturnInterval.py


Input files:
    - ../../data/shapes/cerrado_grid_<grid_id>_cropped.shp
    - ../../data/processed/MBfogo_c<collection>0_cerrado_polygons_b150_CM/MBfogo_c<collection>0_cerrado_polygons_<year>.<csv/shp>
    - ../../data/raw/MBfogo_ba/MBfogo_c<collection>0_cerrado_<year>.tif
    
Output files:
    - ../../data/processed/summary_tables/<nyears>_year_periods/burnedArea_fireCount_intersection_monthFromPixel_grid<grid_id>_period<>.csv
    - ../../data/processed/summary_tables/<nyears>_year_periods/burnedArea_fireCount_intersection_monthFromPolygon_grid<grid_id>_period<>.csv

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
    df_area_cell_monthPolygon = []
    
    # Loop over polygons to assign them to cells individually
    for i in range(len(polygons)): 
        
        # Selecting only the polygon geometry, id and month
        pol_geom = polygons.loc[i, 'geometry']
        pol_month = polygons.loc[i,'month']
        pol_area = polygons_data.loc[i, 'area']
        pol__id = polygons_data.loc[i, 'ID']

        
        if pol_area >= 0.03:
            # First, get the list of cells intersecting the polygon
            cell_inters = grid.loc[grid['geometry'].intersects(pol_geom), 'id'].apply(int).reset_index(drop = True)
            # print(cell_inters)
            
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
                    with open("./exception_50km.txt", "a") as myfile:
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
                    
                    
                    # Now I have the information to create the two dataframes (the only 
                    # difference between dataframes is whether month area comes from pixel 
                    # information or from the month assigned to polygon)
                    
                    # First, the dataframe where month area comes form pixel level information
                    # List of cells that intersect the polygon, along with the total intersection data (all months summed up)
                    area_cell_monthPixel = pd.DataFrame(
                        {'id': cell_inters,
                         'area_T': area_per_geom})
                    # Concatenating the area per month (that comes from pixel-level)
                    area_cell_monthPixel = pd.concat([area_cell_monthPixel, area_per_geom_month], axis = 1)
        
                    
                    # Second, the dataframe where month area comes from the polyon's month
                    area_cell_monthPolygon = pd.DataFrame(
                        {'id': cell_inters,
                         'area_T': area_per_geom})
                    # Generating the month part...
                    # ...creating empty dataframe with all the month columns
                    month_df = pd.DataFrame(columns = [str(m) for m in range(1, 12+1)])
                    # ...assigning the area burned to each cell
                    month_df[str(pol_month)] = area_per_geom
                    # ...filling the rest with 0s
                    month_df = month_df.fillna(0)
                    # ...concatenating to the rest of the dataframe
                    area_cell_monthPolygon = pd.concat([area_cell_monthPolygon, month_df], axis = 1)
        
                        
                    # Appending to list of dataframes
                    df_area_cell_monthPixel.append(area_cell_monthPixel)
                    df_area_cell_monthPolygon.append(area_cell_monthPolygon)
                    
                    # Freeing memory
                    del area_cell_monthPixel, area_cell_monthPolygon, month_df, 
                    area_per_geom, area_per_geom_month, cell_inters
            

                
    # If no polygons were larger than the threshold, or none were intersecting cells
    if len(df_area_cell_monthPixel) > 0:
        # Cast list of dataframes as a unique dataframe for each one of the two options
        df_area_cell_monthPixel = pd.concat(df_area_cell_monthPixel, ignore_index = True)
        df_area_cell_monthPolygon = pd.concat(df_area_cell_monthPolygon, ignore_index = True)
        
        # This dataframe may contain each cell several times (as many times as polygons intersect it)
        # Hence, we can just count the number of times the cell appears as the number of fire polygons, and sum
        # the area of each polygon within cell as total area burned.
        
        # Total area per cell (and month or total)
        temp_monthPixel = df_area_cell_monthPixel.groupby(by = 'id').sum().reset_index()
        temp_monthPolygon = df_area_cell_monthPolygon.groupby(by = 'id').sum().reset_index()
        
        # Counting number of times each cell appears
        df_area_cell_monthPixel = df_area_cell_monthPixel.groupby(by = 'id').size().to_frame('npolygons').reset_index()
        df_area_cell_monthPolygon = df_area_cell_monthPolygon.groupby(by = 'id').size().to_frame('npolygons').reset_index()
        
        # Merging the information
        df_area_cell_monthPixel = df_area_cell_monthPixel.merge(temp_monthPixel, on = 'id', how = 'outer')
        df_area_cell_monthPolygon = df_area_cell_monthPolygon.merge(temp_monthPolygon, on = 'id', how = 'outer')
        
        # Freeing memory
        del temp_monthPixel, temp_monthPolygon
    
    # Otherwise, we'll return empty dataframes
    if len(df_area_cell_monthPixel) == 0:
        df_area_cell_monthPixel = pd.DataFrame(columns = ['id', 'area_T'] + [str(m) for m in range(1, 12+1)])
        df_area_cell_monthPolygon = pd.DataFrame(columns = ['id', 'area_T'] + [str(m) for m in range(1, 12+1)])
        
    
    return df_area_cell_monthPixel, df_area_cell_monthPolygon


"""
MAIN
"""

def main():
    
    # ----------------------- User inputs -------------------------------------
    
    # Grid size
    grid_size_units = '30km'
    
    # UPDATE 24/10/2023
    # MapBiomas fogo collection
    collection = 2
    
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
    
    
    # UPDATE (24/10/2023): 
    # Pixel output filepath
    pixel_output_fp = '../../data/processed/summary_tables/{}_year_periods/burnedArea_fireCount_intersection_monthFromPixel_grid{}_period{}.csv'
    # POlygon output filepath
    polygon_output_fp = '../../data/processed/summary_tables/{}_year_periods/burnedArea_fireCount_intersection_monthFromPolygon_grid{}_period{}.csv'
    
    annual_pixel_output_fp = '../../data/processed/summary_tables/annual/burnedArea_fireCount_intersection_monthFromPixel_grid{}_year{}.csv'
    annual_polygon_output_fp = '../../data/processed/summary_tables/annual/burnedArea_fireCount_intersection_monthFromPolygon_grid{}_year{}.csv'
    
    # ------ Classifying polygon by cell for each year in period -------------- 
    
    
    # Grid filepath
    grid_fp = '../../data/shapes/cerrado_grid_{}_cropped.shp'.format(grid_size_units)
    # Reading grid
    grid = gpd.read_file(grid_fp)
    # Selecting only relevant columns
    grid = grid[['id', '%area_crop', 'withinCerr','geometry']]
    

    
    # Counter variable
    counter = 1
    
    # We work on a group of years at a time
    for (start_year, end_year) in group_years[1:]:
        
        period = '({}, {})'.format(start_year, end_year)
        
        print('Working on group ' + period)
        
        # Creating empty lists where we accumulate the data for all years in group
        df_pixel = []
        df_polygon = []
        
    
        # For each year in group, calculate the number of polygons intersecting 
        # each cell and the total burned area of the polygons (only the area
        # inside the cell)
        for year in range(start_year, end_year + 1):
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
            
            # TODO: erase this later
            # polygon_data = polygon_data[polygon_data['month'].isin([1,2])]
            
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
                output = p.starmap(countCells_area_intersectingPolygon, pool_input)
                df_pixel_year = [r[0] for r in output]
                df_polygon_year = [r[1] for r in output]
                
                # Freeing memory
                del output

                print('Finished parallel computation.')
                
                
            # Converting list of dataframes to a single dataframe
            df_pixel_year = pd.concat(df_pixel_year, ignore_index = True)
            df_polygon_year = pd.concat(df_polygon_year, ignore_index = True)
            
            # These dataframes contain the total area and area per month for each cell
            # But cells may be duplicated in the dataframe (polygons in different slices intersected the same cell)
            # So, we get sum to get the total count and area
            df_pixel_year = df_pixel_year.groupby(by = 'id').sum().reset_index()
            df_polygon_year = df_polygon_year.groupby(by = 'id').sum().reset_index()
            
            # -----------------------------------------------------------------
            # UPDATE 17/11/2023: I need the annual information per cell!!
            
            # Create column with year id
            df_polygon_year['year'] = year
            # Casting columns in order
            df_polygon_year = df_polygon_year[['id', 'year','npolygons','area_T'] + [str(m) for m in range(1, 12+1)]]
            # Saving to file
            df_polygon_year.to_csv(annual_polygon_output_fp.format(grid_size_units, year), index = False)
            # Dropping year column
            df_polygon_year = df_polygon_year.drop(columns = 'year')
            
            
            # Create column with year id
            df_pixel_year['year'] = year
            # Casting columns in order
            df_pixel_year = df_pixel_year[['id', 'year','npolygons','area_T'] + [str(m) for m in range(1, 12+1)]]
            # Saving to file
            df_pixel_year.to_csv(annual_pixel_output_fp.format(grid_size_units, year), index = False)
            # Dropping year column
            df_pixel_year = df_pixel_year.drop(columns = 'year')
            
            
            # -----------------------------------------------------------------
            
            # Now we have all the data for each cell in year, we accumulate it to next year
            df_pixel.append(df_pixel_year)
            df_polygon.append(df_polygon_year)
            
        
            
            # Freeing memory
            del df_pixel_year, df_polygon_year
            
        print('Formatting data of all years in period...')
            
        # Now we have the data for all years in group as two lists of dataframes
        # First we convert them to dataframes
        df_pixel = pd.concat(df_pixel, ignore_index = True)
        df_polygon = pd.concat(df_polygon, ignore_index = True)
        
        # Finally, again, these dataframes contain duplicated cells (one time per each
        # year in group where they have had a fire). Summarise
        df_pixel = df_pixel.groupby(by = 'id').sum().reset_index()
        df_polygon = df_polygon.groupby(by = 'id').sum().reset_index()
        
        
        # Now we add a column to each indicating the period
        df_pixel['period'] = period
        df_polygon['period'] = period
        
        # Reshuffling columns
        df_pixel = df_pixel[['id','period','npolygons','area_T'] + [str(m) for m in range(1, 12+1)]]
        df_polygon = df_polygon[['id','period','npolygons','area_T'] + [str(m) for m in range(1, 12+1)]]
        
        print('Saving data for all years in period...\n')
        
        # Saving dataframes
        df_pixel.to_csv(pixel_output_fp.format(n_years, grid_size_units, counter), index = False)
        df_polygon.to_csv(polygon_output_fp.format(n_years, grid_size_units, counter), index = False)
        
        counter += 1
        
        print('--------------------------------------------------------------')
        print('\n\n')
        
        # brea!k
    

        
        
if __name__ == '__main__':
    main()

