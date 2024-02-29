# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 10:25:41 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculate the median and 95th quantile of each topographic measure
(roughness, Topographic Ruggedness Index, and slope) per cell. It carries out
these calculations for all cells in grid, regardless of their propoertion of 
area inseide the Cerrado.

"""

"""
Imports
"""
import numpy as np
import rasterio
from rasterio.mask import mask
import geopandas as gpd
import pandas as pd

from pytictoc import TicToc
t = TicToc()

from tqdm import tqdm



"""
Functions
"""
def get_grid_rasterCRS(grid_fp, raster_fp):
    
    # Get the raster's CRS
    with rasterio.open(raster_fp) as src:
        raster_crs = src.crs
        
    
    # Read the grid
    grid = gpd.read_file(grid_fp)
    
    # Transform to rsater's CRS
    grid = grid.to_crs(raster_crs)
    
    return(grid)

def topographicMeasure_summary(raster_fp, polygon):
    
    # Open raster in read mode
    with rasterio.open(raster_fp) as raster:
        
        # Read only the raster information within the polygon
        out_img, out_transform = mask(dataset = raster, shapes=[polygon], crop=True, filled = False)
        
    # Getting non-masked values
    out_img = out_img.compressed()
    
    # Calculating median and 95th percentile
    topographic_median = np.percentile(out_img, 50)
    topographic_q95 = np.percentile(out_img, 95)

    
    return(topographic_median, topographic_q95)


"""
MAIN
"""

def main():
    
    # --------------------- User inputs --------------------------------------
    
    grid_id = '30km'
    
    # Topographic measures to work with
    topo_measures = ['roughness', 'TRI', 'slope']
    
    # Clipped grid (so that we only consider topographic measures within the 
    # Cerrado)
    grid_fp = '../../data/shapes/cerrado_grid_{}_cropped_partlyClipped.shp'.format(grid_id)
    
    # Raster file with the topographic measures
    topo_raster_fp = '../../data/processed/SRTM_DEM/SRTM_DEM_merged_{}.tif'
    
    # Raster file with elevation topographic measure
    elev_raster_fp = '../../data/raw/SRTM_DEM/SRTM_DEM_merged.tif'
    
    # Output file to save the calculations
    output_fp = '../../data_new_grids/processed/summary_tables/topography_allVars_grid{}.csv'.format(grid_id)
    
    
    # --------------------- Processing the data ------------------------------
    
    grid = get_grid_rasterCRS(grid_fp, topo_raster_fp.format(topo_measures[0]))
    
    
    # Create list of dataframes to store each cell's topographic data
    topo_summary = []
    
    # Working on one cell at a time
    t.tic()
    for cell in tqdm(grid['id']):
        
        
        # Create dataframe to store the cells' topographic data
        cell_data = pd.DataFrame({'id': [cell]})
        
        # Fetching the cell's geometry
        cell_geom = grid.loc[grid['id'] == cell, 'geometry'].item()
        
        # Calculating the summary statistics for each topographic measure at
        # a time
        for measure in topo_measures:
        
            # Input filepath for measure
            raster_fp = topo_raster_fp.format(measure)
            
            cell_data[measure + '_q50'], cell_data[measure + '_q95'] = \
                topographicMeasure_summary(raster_fp, cell_geom)
                
        # UPDATE (06/10/2023): adding elevation topographic measure
        cell_data['elevation_q50'], cell_data['elevation_q95'] = \
            topographicMeasure_summary(elev_raster_fp, cell_geom)
       
        
        # Append cell's data to list of dataframes
        topo_summary.append(cell_data)
    
    t.toc("Cell's topographic data calculated in")
        
    
    # Convert list of dataframes into unique dataframe
    topo_summary = pd.concat(topo_summary, ignore_index = True)
    
    # Save data
    print('\nSaving data...')
    topo_summary.to_csv(output_fp)
    

    
        
        
        
        
        

    

if __name__ == '__main__':
    main()

