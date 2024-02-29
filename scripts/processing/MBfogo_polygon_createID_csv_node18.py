# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 13:17:18 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that wrangles a GeoDataFrame into the format:
    <year><month>_<dataframe_id> | year | month | geometry
and overwrites it. It also generates the same information - without the geometry 
column - and saves it to a CSV.

Input files:
    ../../data/processed/MBfogo_c10_cerrado_polygons_b150_CM/MBfogo_c10_cerrado_polygons_<year>.shp
    
Output files:
    ../../data/processed/MBfogo_c10_cerrado_polygons_b150_CM/MBfogo_c10_cerrado_polygons_<year>.shp
    ../../data/processed/MBfogo_c10_cerrado_polygons_b150_CM/MBfogo_c10_cerrado_polygons_<year>.csv    
    
"""

"""
Imports
"""
import geopandas as gpd

from multiprocessing import Pool

import tqdm
import sys
sys.path.append('../utils/') 
import fileList

from MBfogo_extract_time_reference import extract_year
from create_yearMonthID import create_yearMonthID

"""
Functions
"""

def cleanData(polygon_fp):
    """
    Function that receives a filepath to a GeoDataFrame, reads it, gives an ID 
    to each polygon and formats it into:
        ID | year | month | geometry
    It also creates a CSV file with the same information (except the geometry).

    Parameters
    ----------
    polygon_fp : string
        Filepath to the GeoDataFrame that neds to be formatted.

    Returns
    -------
    None. Outputs are SHP and CSV files that are saved at the end of this function.

    """
    
    # 1. Read GeoDataFrame
    polygon_gdf = gpd.read_file(polygon_fp)
    
    # UPDATE (06/10/2023)
    # Renaming column raster_val to month
    polygon_gdf = polygon_gdf.rename(columns = {'raster_val':'month'})

    
    # 2. Create year columns
    # Extract year and month from filename
    year = extract_year(polygon_fp)
    polygon_gdf['year'] = int(year)
    # If polygon already had column year, this will be overwritten. Hence, it 
    # is important that the polygon_fp ends in _year.shp
    
    # UPDATE 30-03-2023: this is not necessary because now the month is assigned
    # by ./scripts/processing/MBfogo_assignMonth_toPolygon.py
    # 3. Rename column 'raster_val' to 'month'
    # polygon_gdf = polygon_gdf.rename(columns = {'raster_val': 'month'})
    # This line will only be effective if there is a column raster_val in GeoDataFrame.
    # Otherwise, it will hav eno effect.
    
    
    # 4. Reordering
    polygon_gdf = polygon_gdf[['year', 'month','geometry']]
    
    
    # 5. Create unique ID as column 
    polygon_gdf = create_yearMonthID(polygon_gdf)
    
    # The ID is returned as the index, cast into column
    polygon_gdf = polygon_gdf.reset_index()

    
    # 6. Save index column, year and month as a CSV file
    polygon_gdf.loc[:, polygon_gdf.columns != 'geometry'].to_csv(polygon_fp[:-4] + '.csv', index=False)
    
    # Save reformatted dataframe overwriting the previous one
    polygon_gdf.to_file(polygon_fp)
    



"""
MAIN
"""

def main():
    
    # ------------------ User inputs -----------------------------------------
    
    # Read list of filepaths from directory
    # in_dir = '../../data/processed/MBfogo_c10_cerrado_polygons_b150_CM/'
    in_dir = '../../data/processed/MBfogo_c20_cerrado_polygons/'
    list_fp =fileList.fileList(in_dir, in_extensions = '.shp')
    # list_fp = [in_dir + 'MBfogo_c10_cerrado_polygons_1986.shp',
    #            in_dir + 'MBfogo_c10_cerrado_polygons_1987.shp',
    #            in_dir + 'MBfogo_c10_cerrado_polygons_1988.shp',
    #            in_dir + 'MBfogo_c10_cerrado_polygons_1989.shp']   
    list_fp = fileList.selecting_fp_withCondition(list_fp, [str(i) for i in range(2015, 2022 + 1)])
    
    # This program can be run sequentially over list_fp or parallelising
    # over the list of filepaths
    # option = 'sequential'
    option= 'parallel'
    
    # number of processes to use if running in parallel
    num_proc = 5
    
    # --------------- Formatting GeoDataFrame --------------------------------
    
    # Sequentially
    if option == 'sequential':
        for fp in list_fp:
            print('Working on file ' + fp)
            cleanData(fp)
            
    # In parallel
    elif option == 'parallel':
        with Pool(num_proc) as p:
            print('Working in parallel...')
            tqdm.tqdm(p.map(cleanData, list_fp))
    
if __name__ == '__main__':
    main()


