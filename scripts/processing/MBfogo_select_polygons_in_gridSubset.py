# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 11:54:39 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Brief script that reads a subset of cells of a grid (a shapefile and a csv per year),
extracts the list of cells, and uses it to select MapBiomas FOGO polygons from the files.
"""

"""
Imports
"""
import geopandas as gpd
import pandas as pd
from pytictoc import TicToc
t = TicToc()

"""
MAIN
"""

if __name__ == '__main__':
    
    # ----------------- User inputs -----------------------------
    
    # list of years to work on 
    list_years = [1990, 2000, 2010, 2020]
    
    # Input file
    # in_dir = '../../data/processed/MBfogo_c10_cerrado_polygons/'
    # in_file = 'MBfogo_c10_cerrado_polygons_{}.{}'
    in_dir = '../../data/tests/fire_polygons/nobuff_consecutiveMonths/'
    in_file = 'MBfogo_c10_cerrado_polygons_subset_{}_nobuff.{}'
    
    # Output directory
    # out_dir = '../../data/tests/fire_polygons/individual/'
    out_dir = in_dir
    
    # The subset of cells
    grid_subset_fp = '../../data/tests/fire_polygons/test_area.shp'
    
    # ----------------- Subsetting cells ------------------------
    
    # Reading list of cells
    grid_subset = gpd.read_file(grid_subset_fp)
    # List of cells
    grid_subset = grid_subset['id']
    
    # Subset the shapefile and the csv file on a year by year basis
    for year in list_years: 
        
        print('Working on year {}'.format(year))
        
        # First, read both files
        csv_polygons = pd.read_csv(in_dir + in_file.format(year, 'csv'))
        
        # Subset polygons in csv if assigend to 200 km cell in area subset
        csv_polygons = csv_polygons[csv_polygons['cellID_200km'].isin(grid_subset)]
        
        
        # Then, read polygons and subset those in csv list
        print('Reading SHP file...')
        t.tic()
        shp_polygons = gpd.read_file(in_dir + in_file.format(year, 'shp'))
        t.toc()
        # Subset
        shp_polygons = shp_polygons[shp_polygons['ID'].isin(csv_polygons['ID'])]
        
        # Write both files to ouput dir
        csv_polygons.to_csv(out_dir + in_file.format(year, 'csv'))
        print('Writing SHP file...')
        shp_polygons.to_file(out_dir + in_file.format(year, 'shp'))

        # break
        print('')
        
    

