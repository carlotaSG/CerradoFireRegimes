# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 15:08:50 2022

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculates the distance matrix for all cell pairs in a grid. It
retrieves the distance matrix in condensed form.

Input data:
    - ../../data/shapes/cerrado_grid_<grid_id>_cropped.shp
    
Output data:
    - ../../data/processed/clustering_algorithm/distance_matrix_grid_<grid_id>_epsg5880.csv

Written with inspiration from:
    - (08/08/2022): https://stackoverflow.com/questions/36696613/calculating-distance-between-multiple-sets-of-geo-coordinates-in-python
    - (08/08/2022): https://www.faqcode4u.com/faq/170881/melt-the-upper-triangular-matrix-of-a-pandas-dataframe
    

"""

"""
Imports
"""
import geopandas as gpd
import geopy.distance
import numpy as np
import pandas as pd
from pytictoc import TicToc
t = TicToc()
from itertools import combinations_with_replacement
from multiprocessing import Pool


"""
Functions
"""
def calculate_distanceMatrix(df):
    # Function that creates distance matrix from list of coordinates
    
    # Create column with coordinates as tuple
    df['latlon'] = list(zip(df['lat'], df['lon']))
    
    # Create empty squared dataframe to store the distances between each pair of
    # coordinates (one column and row per coordinate)
    square = pd.DataFrame(
        np.zeros((df.shape[0], df.shape[0])),
        index=df.index, columns=df.index
    )
    


    # Function within function that calculates the distances between each pair
    # coordinates
    def get_distance(col):
        
        end = df.loc[col.name, 'latlon']
        distance =  df['latlon'].apply(geopy.distance.distance,
                                  args=(end,),
                                  ellipsoid='WGS-84'
                                 )
        distance = distance.apply(lambda x: x.km)
        return distance
    
    
    # Calculate distances between cells (pairs of coordinates)
    distances = square.apply(get_distance, axis=1).T
    
    # Change row and column names to those of the cell's id
    grid['id2'] = grid['id'] # Creating dummy column to avoid later confusion
    distances.index = grid['id']
    distances.columns = grid['id2']
    
    return distances
    
    
def longFormat(df):
    
    # Turn to NaNs the lower triangle
    df = df.where(np.triu(np.ones(df.shape)).astype(bool))
    
    
    # Cast in long format
    df = df.stack()
    # print(df)
    # df = df.rename(columns = {'id':'id2'}, level = 1)
    df = df.reset_index()

    # Format dataframe
    df.columns = ['id1','id2','distance_km']
    
    return df

def calculateDistance(gdf1, gdf2):
    
    gdf1['dist_m'] = gdf1.distance(gdf2['geometry'], align = False)
    
    # Output dataframe
    df = gdf1[['id1','dist_m']].join(gdf2[['id2']])
    
    return df

"""
MAIN
"""

if __name__ == '__main__':
    
    # User inputs ------------------------------------------------------------
    
    # Grid identifier
    grid_id = '50km'
    
    # # Theshold to impose to cells: we only work with cells whose area in the 
    # # Cerrado is larger than x%
    # cell_thresh = 75.0
    
    # # Extra cells to include and cells to exclude
    # include_cells = [567, 568]
    # exclude_cells = [1451]
    
    
    # Input grid shapefile
    # grid_fp = '../../data/shapes/test_80cells_{}.shp'.format(grid_id)
    grid_fp = '../../data/shapes/cerrado_grid_{}_cropped.shp'.format(grid_id)
    
    # Output file
    out_fp = '../../data/processed/clustering_algorithm/distance_matrix_grid_{}_epsg5880.csv'.format(grid_id)
    # out_fp = '../../data/processed/clustering_algorithm/distance_matrix_80cells_grid_{}_epsg5880.csv'.format(grid_id)
    
    
    # Read data and format ---------------------------------------------------

    
    grid = gpd.read_file(grid_fp)
    # grid = grid[((grid['%area_crop'] >= cell_thresh) | grid['id'].isin(include_cells)) & ~grid['id'].isin(exclude_cells)].reset_index(drop = True)
    
    crs = grid.crs
    
    # Selecting only relevant columns
    grid = grid[['id','lat','lon']]
    
    # Then getting the list of cells
    cells = grid['id'].unique().tolist()
    
    # Generating all pairs of cells possible in the data
    cells = pd.DataFrame(list(combinations_with_replacement(cells, 2)), columns = ['id1','id2'])
    
    # From this dataframe, create two geodataframes, each with their own geometries (points from lat lon)
    gdf1 = cells[['id1']].merge(grid, how = 'left', left_on = 'id1', right_on = 'id').drop(columns = 'id')
    gdf1 = gpd.GeoDataFrame(gdf1, geometry = gpd.points_from_xy(gdf1.lon, gdf1.lat), crs = crs)
    
    # Making sure the CRS is projected
    gdf1 = gdf1.to_crs('epsg:5880')
    
    gdf2 = cells[['id2']].merge(grid, how = 'left', left_on = 'id2', right_on = 'id').drop(columns = 'id')
    gdf2 = gpd.GeoDataFrame(gdf2, geometry = gpd.points_from_xy(gdf2.lon, gdf2.lat), crs = crs)
    
    # Making sure the CRS is projected
    gdf2 = gdf2.to_crs('epsg:5880')
    
    # Freeing memory
    del grid
    
    
    # Calculate distance -----------------------------------------------------
    
    # Preparing data to calculate distances in parallel
    
    # Number of chunks
    nchunks = 28
    # Number of elements per chunk
    nelements = round(len(gdf1) / nchunks)
    
    input_toParallel = [(gdf1.iloc[i : i+nelements, :].reset_index(drop = True), gdf2.iloc[i : i+nelements, :].reset_index(drop = True)) for i in range(0, len(gdf1), nelements)]
    
    
    # Freeing memory
    del gdf1, gdf2
    
    t.tic()
    with Pool(nchunks) as p: 
            
            print("Starting parallel computation in calculateDistance...")
            
            # output = p.starmap(assign_polygon_to_cell, pool_input)
            out_df = p.starmap(calculateDistance, input_toParallel)
    t.toc('Parallisation finished in')
            
    # Casting list of dataframes as unique dataframe
    out_df = pd.concat(out_df, ignore_index = True)
    
    # Reordering columns
    out_df = out_df[['id1','id2','dist_m']]
    
    out_df.to_csv(out_fp, index = False)
    