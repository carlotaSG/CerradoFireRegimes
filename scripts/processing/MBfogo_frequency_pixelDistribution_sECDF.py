# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 15:23:23 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that gets the distribution of burned frequency at pixel level
"""

"""
Imports
"""
import geopandas as gpd
import pandas as pd
import numpy as np
import numpy.ma as ma
import rasterio
from rasterio.mask import mask
# from statsmodels.distributions.empirical_distribution import ECDF
import statsmodels.distributions.empirical_distribution as edf
from scipy.interpolate import interp1d

from pytictoc import TicToc
t = TicToc()
from multiprocessing import Pool

"""
Functions
"""
def getFrequencies_fromPixel(cell_geometry, cell_id, start_year, end_year,
                             mbfogo_fp, mblulc_fp, lulc_col):
    
    
    # List to accumulate the rasters with the burned/unburned and flammable pixels
    # pixels for each different year
    burned_flammable = []
    for year in range(start_year, end_year + 1):
        
        # print(year)
        
        # Read the burned area raster using the cell geometry to mask the raster
        # Open raster file in 'r' mode
        with rasterio.open(mbfogo_fp.format(year)) as raster:
            # print('here2')
            mbfogo, _ = mask(dataset = raster, shapes = [cell_geometry], crop = True, filled = False)
            
            # removing empty dimensions
            mbfogo = np.squeeze(mbfogo)
            
            # Converting burned pixels to 1
            temp = np.where(mbfogo >=1 , 1, 0)
            
            # Returning masked pixels to masked
            mbfogo = ma.masked_where(np.isnan(mbfogo), temp)
            mbfogo = np.where(mbfogo.mask, np.nan, mbfogo)
            
            
            # Freeing memory
            del temp
            
            
        
        # Now we will cast all non-flammable pixels as nans to take them out of
        # the frequency calculation
        # Read the land-use data raster using the cell geometry to mas the raster
        with rasterio.open(mblulc_fp.format(year)) as raster:
            
            mblulc, _ = mask(dataset = raster, shapes = [cell_geometry], crop = True, filled = False)
            
            # removing empty dimensions
            mblulc = np.squeeze(mblulc)

            # UPDATE (24/10/2024): including the possibility of changing the MapBimas lulc collection
            # Converting non-flammable pixels to nans
            if lulc_col == 6:
                # nonflammable (comprises: Salt Flat (32), Rocky Outcrop (29), 
                #                 Non vegetated area (22), Water(26), Non Observed (27)) -> 'nonflammable'
                temp = np.where(np.isin(mblulc, [32, 29, 23, 24, 30, 25, 33, 31, 27]),
                                np.nan,
                                1)
            elif lulc_col == 8:
                # nonflammable (comprises: Salt Flat (32), Rocky Outcrop (29), 
                #                 Non vegetated area (22), Water(26), Non Observed (27)) -> 'nonflammable'
                temp = np.where(np.isin(mblulc, [32, 29, 50, 23, 24, 30, 25, 33, 31, 27]),
                                np.nan,
                                1)
                
            # Returning masked pixels to masked
            mblulc = ma.masked_where(np.isnan(mblulc), temp)
            # and converting masked to nans
            mblulc = np.where(mblulc.mask, np.nan, mblulc)
            
            # Returning masked pixels to masked
            mblulc = ma.masked_where(np.isnan(mblulc), temp)
            # and converting masked to nans
            mblulc = np.where(mblulc.mask, np.nan, mblulc)
            

            
        # Multiply the flammable pixels and the mbfogo pixels to discard those 
        # burned and unburned pixels that are non-flammable in this year
        burned_flammable_year = mbfogo * mblulc

        
        burned_flammable.append(burned_flammable_year)
        
    
    # Finally, for each year, we know if a pixel was non-flammable at some point
    # (it will be ignored in further calculations) and wether it burned in any year
    # Let's calculate how many times each pixel that has never been non-flammable 
    # burned in the period
    burned_flammable = sum(burned_flammable)
    
    
    # Flatten the array 
    burned_flammable = burned_flammable.flatten()
    
    # Discard nans
    burned_flammable = burned_flammable[~np.isnan(burned_flammable)]

    
    return(burned_flammable)


def calculate_sECDF(data, n):
    
    # First, fit the ecdf
    data_ecdf = edf.ECDF(data)
    

    slope_changes = list(range(0, n + 1))
    
    data_ecdf_values_at_slope_changes = [data_ecdf(item) for item in slope_changes]
    
    # print(data_ecdf_values_at_slope_changes)
    
    
    
    inverted_ecdf = interp1d(data_ecdf_values_at_slope_changes,
                             slope_changes)
    
    # It may be that the 5th quantile is below the interpolation range's minimum
    # value. In that case, we use this value instead of the 5th quantile
    if 0.05 < data_ecdf_values_at_slope_changes[0]:
        q05val = data_ecdf_values_at_slope_changes[0]
    else:
        q05val = 0.05
        
    if 0.5 < data_ecdf_values_at_slope_changes[0]:
        q50val = data_ecdf_values_at_slope_changes[0]
    else:
        q50val = 0.5
        
    if 0.95 < data_ecdf_values_at_slope_changes[0]:
        q95val = data_ecdf_values_at_slope_changes[0]
    else:
        q95val = 0.95
        
    if 0.99 < data_ecdf_values_at_slope_changes[0]:
        q99val = data_ecdf_values_at_slope_changes[0]
    else:
        q99val = 0.99
    
    q05 = inverted_ecdf(q05val)
    q50 = inverted_ecdf(q50val)
    q95 = inverted_ecdf(q95val)
    q99 = inverted_ecdf(q99val)
    
    return q05, q50, q95, q99


def getFrequencies_inParallel(cell_geometry, cell_id, period,
                              mbfogo_fp, mblulc_fp, nyears_period, lulc_col):
    
    # Getting the array of frequencies for each non-flammable pixel in the period
    frequencies = getFrequencies_fromPixel(cell_geometry, cell_id, period[0], period[1],
                                           mbfogo_fp, mblulc_fp, lulc_col)
    
    
    # Cast this information as a dataframe
    df_freq_cell_period = pd.DataFrame({'id': cell_id, 'period': str(period), 'frequency': frequencies})
    
    # Summarise
    df_freq_cell_period = df_freq_cell_period.groupby(by = ['id','period','frequency']).size().reset_index(name = 'npixels')
    
    # If there is only one frequency (for some reason), then the quantiles are just this frequency
    if len(df_freq_cell_period) == 1:
        
        the_freq = df_freq_cell_period.loc[0, 'frequency'].item()
        f05, f50, f95, f99 = the_freq, the_freq, the_freq, the_freq
        
    else:
    
        # Calculate the quantiles of the fire frequency
        f05, f50, f95, f99 = calculate_sECDF(frequencies, nyears_period)
    
    # # print(df_freq_cell_period)
    
    # # Let's get the quantiles of the frequency distribution
    # f05, f50, f95, f99 = getQuantiles(frequencies, nyears_period)
    
    # Casting as dataframe
    df_quant_cell_period = pd.DataFrame({'id': [cell_id], 'period': str(period), 
                             'f05': f05, 'f50': f50, 
                             'f95': f95, 'f99': f99})
    
    return df_freq_cell_period, df_quant_cell_period
    

"""
MAIN
"""
def main():
    
    # ----------------------- User inputs -------------------------------------
    
    grid_id = '30km'
    
    grid_fp = '../../data/shapes/cerrado_grid_{}_cropped.shp'.format(grid_id)
    
    periods = [(1985, 1993), (1994, 2002), (2003, 2011), (2012, 2020)]
    
    nperiods = len(periods)
    
    nyears_period = 9
    
    # UPDATE (24/10/2023): Including all cells
    # # Threshold area
    # area_thresh = 75.0
    
    # # Extra cells to include and cells to exclude
    # include_cells = [567, 568]
    # exclude_cells = [1451]
        
    # UPDATE (24/10/2023):  Updating the MapBiomas land-use and fogo collection
    col_fogo = 2
    col_lulc = 8
    
    mbfogo_fp = '../../data/raw/MBfogo_ba/MBfogo_c{}0_cerrado_{}.tif'.format(col_fogo, '{}')
    
    # mblulc_fp = 'D:/carlota/projects/data/MBlanduse_c06/mapbiomas-brazil-collection-60-cerrado-{}.tif'
    # mblulc_fp = '../../data/raw/MBlulc/mapbiomas-brazil-collection-60-cerrado-{}.tif'
    mblulc_fp = '../../data/raw/MBlulc/mapbiomas-brazil-collection-80-{}.tif'
    
    num_process = 38
    
    
    freq_fp = '../../data/processed/summary_tables/{}_year_periods/fireFrequency_fromPixels_sECDF_grid{}_period{}.csv'
    
    quant_fp = '../../data/processed/summary_tables/{}_year_periods/fireFrequency_fromPixels_sECDF_quantiles_grid{}_period{}.csv'
    
    
    # ---------------------- Reading grid -------------------------------------
    
    grid = gpd.read_file(grid_fp)
    
    # UPDATE (24/10/2023): Including all cells
    # # Selecting cells with a certain area inside the Cerrado
    # grid = grid[((grid['%area_crop'] >= area_thresh) | grid['id'].isin(include_cells)) & ~grid['id'].isin(exclude_cells)].reset_index(drop = True)
    
    # Convert the grid's CRS to that of MapBiomas'
    grid = grid.to_crs('epsg:4326')
    
    
    # ------------------ Calculating frequencies per period -------------------
    
    
    for i, period in zip(range(1, nperiods + 1), periods):
        
        print('Working on period {}...'.format(period))
        
        # df_freq = []
        # df_quant = []
        
        # This will be parallelised
        
        # Preparing input for parallel computation
        
        # Number of cells
        ncells = len(grid)
        t.tic()
        print('Starting parallel computation...')
        poolInput = list(zip([features for features in grid['geometry']],
                             grid['id'],                             
                             [period]*ncells, 
                             [mbfogo_fp]*ncells,
                             [mblulc_fp]*ncells,
                             [nyears_period]*ncells,
                             [col_lulc]*ncells))
        
        
        with Pool(num_process) as p:
            # Sending to compute in parallel
            output = p.starmap(getFrequencies_inParallel, poolInput)

        df_freq_period = [output[i][0] for i in range(len(output))]
        df_quant_period = [output[i][1] for i in range(len(output))]
        
        df_freq_period = pd.concat(df_freq_period, ignore_index = True)
        df_quant_period = pd.concat(df_quant_period, ignore_index = True)

        
        print('Casting together all period data...')

        
        df_freq_period.to_csv(freq_fp.format(nyears_period, grid_id, i), index = False)
        df_quant_period.to_csv(quant_fp.format(nyears_period, grid_id, i), index = False)
        
        t.toc('Period {} done in'.format(period))

    
if __name__ == '__main__':
    main()