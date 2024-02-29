# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 17:15:51 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculates the average value per period for various climatic variables.

For each year it reads the average value for each variable, cell, year and month.
Then it calculates the year's average, and then the period's average.

"""

"""
Imports
"""
import numpy as np
import pandas as pd


import sys
sys.path.append('../utils/')
from readFormat_dataYears import readFormat_dataYears
from tqdm import tqdm

"""
Functions
"""
def MCWD(data, year):
    # Function that gets the CWD for a 12-month period and returns the MCWD
    # along with the corresponding month and year
    
    # Identify the MCWD
    # print(data['CWD'].idxmin())
    mcwd = data.loc[[data['CWD_aet'].idxmin()], :]
    # print(mcwd)
    
    # Rename CWD column to MCWD
    mcwd = mcwd.rename(columns = {'CWD_aet': 'MCWD_aet'})
    
    # Adding column with climate year identifier
    mcwd['climate_year'] = year
    
    return mcwd
    
    
def MCWD_timeseries(data, w_month):
    # Function that gets a time series of CWD and information on the wettest month
    # on average and gets the MCWD for every 12-month period
    
    # Get index positions of wettest months
    w_index = data[data['month'] == w_month].index
    
    # Then, we loop over every 12-month period from the wettest month finding
    # the corresponding MCWD and the corresponding month-year 
    # List to accumulate the MCWD per climate year
    mcwd_data = []
    
    # Initialise counter of climate years
    climate_year = 0
    for i in w_index:     
        
        # List of indices (rows) to work on
        subset_index = list(range(i, i + 12))
        # print(subset_index)
        
                
        # Working with subset of 12 months only if all of them are within the 
        # time series
        if any(x >= (len(data) -1) for x in subset_index):
            break
    
        else:
        
            # print(data.iloc[subset_index, :])
            
            # Subset of 12 months for cell
            data_subset = data.iloc[subset_index, :]
            
            # Calculate MCWD for a subset of 12 months
            mcwd = MCWD(data_subset, climate_year)
            # print(mcwd)
            # Calculate length of the dry season in terms of the MCWD
            mcwd['CWD_length'] = len(data_subset[(data_subset['mean_prec'] - data_subset['mean_aet']) < 0])
            
            # Appending dataframe to list
            mcwd_data.append(mcwd)
            
            
            
            climate_year += 1
    
    mcwd_data = pd.concat(mcwd_data, ignore_index = False)
    
    # Reordering columns
    mcwd_data = mcwd_data[['id', 'year', 'month', 'climate_year', 'MCWD_aet', 'CWD_length']]
    
    return mcwd_data

def identify_wettestmonth(data):
    # Function data has id  | period | year | month | mean_prec | ...
    
    # Subset relevant information
    data = data[['id','period','month','mean_prec']]
    
    # Calculate average precipitation per month
    data = data.groupby(by = ['id','period','month']).mean().reset_index()
    
    # print(data)
    
    # The wettest month is that with the largest precipitation
    wmonth = data.loc[data['mean_prec'].idxmax(), 'month'].item()
    
    return wmonth
    

def CWD(data, et_var):
    # Function that receives a list of 12 month data on precipitation and aet
    # and calcualtes the CWD for each month
    
    # First, forcing the indices of the 12-month data to be from 0 to 11
    data = data.reset_index(drop = True)
    
    for i in range(1, 12):
        # ERA 5 Land provides ET with negative sign, whilst Terraclimate provides it in positive
        # So we take the absolute value of ET
        deficit = data.loc[i-1, 'CWD_' + et_var] + data.loc[i, 'mean_prec'] - abs(data.loc[i, 'mean_' + et_var])
        
        if deficit < 0:
            data.loc[i, 'CWD_' + et_var] = deficit
        elif deficit >= 0:
            data.loc[i, 'CWD_' + et_var] = 0
        
    return data



def CWD_timeseries(data, w_month, et_var):
    # Function that calculates the CWD for a time series of precipitation and 
    # aet per month and year given the wettest month 
    
    
    # First, create column to store the CWD filled with NAs
    data['CWD_' + et_var] = np.nan
    # print(data)
    
    # Then, set CWD for w_month
    data.loc[data['month'] == w_month, 'CWD_' + et_var] = 0.0
    # print(data)
    # Get index positions of wettest months
    w_index = data[data['month'] == w_month].index
    # print(w_index)
    # We loop over sets of 12 months, calculating the CWD with respect to the 
    # wettest month in the subset
    data_CWD = []
    for i in w_index:     
        
        # List of indices (rows) to work on
        subset_index = list(range(i, i + 12))
        # print(subset_index)
        
                
        # Working with subset of 12 months only if all of them are within the 
        # time series
        if any(x >= (len(data) -1) for x in subset_index):
            break
    
        else:
        
            # print(data.iloc[subset_index, :])
            
            # Calculate CWD for a subset of 12 
            data_CWD.append(CWD(data.iloc[subset_index, :], et_var))
            # break
    
    data_CWD = pd.concat(data_CWD, ignore_index = False)
    
    # Merging back to data, to recover the full time series
    # First, keep only CWD column
    data_CWD = data_CWD[['id', 'year', 'month', 'CWD_' + et_var]]
    data_CWD = data.drop(columns = 'CWD_' + et_var).merge(data_CWD, on =['id','year','month'], how = 'left')
    # print(data_CWD)
    return data_CWD


"""
MAIN
"""
def main():
    
    # ------------------------- User inputs ---------------------------------
    
    # Data files and the columns we want from every file
    in_dir = '../../data/processed/climate_data/'
    
    chirps_fp = in_dir + 'chirps_grid50km_{}.csv'
    chirps_cols = ['mean_prec']
    
    era5_fp = in_dir + 'era5land_grid50km_{}.csv'
    era5_cols = ['mean_temp', 'mean_rh', 'mean_vpd']
    
    terraclim_fp = in_dir + 'terraclimate_grid50km_{}.csv'
    terraclim_cols = ['mean_pet', 'mean_aet']
    
    # Year period
    start_year = 1985
    end_year = 2020
    
    # Grid to use
    grid_id = '50km'
    
    # Number of periods we are working with
    n_periods = 4
    
    # Threshold area
    area_thresh = 75.0
    
    # Extra cells to include and cells to exclude
    include_cells = [567, 568]
    exclude_cells = [1451]
    
    output_fp_yearly = in_dir + 'summary_yearly_allVariables_{}periods.csv'.format(n_periods)
    output_fp_period = in_dir + 'summary_period_allVariables_{}periods.csv'.format(n_periods)
    
    # -------------------- Reading and formatting data -----------------------
    
    # Reading and formatting the data
    chirps = readFormat_dataYears(chirps_fp, start_year, end_year, grid_id, chirps_cols,
                             periods = True, nperiods = n_periods, 
                             drop_geometry = True,
                             thresh = area_thresh, incl_cells = include_cells, excl_cells = exclude_cells)
    chirps = chirps.drop(columns = '%area_crop').reset_index(drop = True)
    # print(chirps)
    
    era5 = readFormat_dataYears(era5_fp, start_year, end_year, grid_id, era5_cols,
                             periods = True, nperiods = n_periods, 
                             drop_geometry = True,
                             thresh = area_thresh, incl_cells = include_cells, excl_cells = exclude_cells)
    era5 = era5.drop(columns = '%area_crop').reset_index(drop = True)
    # print(era5)
    
    terraclim = readFormat_dataYears(terraclim_fp, start_year, end_year, grid_id, terraclim_cols,
                             periods = True, nperiods = n_periods, 
                             drop_geometry = True,
                             thresh = area_thresh, incl_cells = include_cells, excl_cells = exclude_cells)
    terraclim = terraclim.drop(columns = '%area_crop').reset_index(drop = True)
    # print(terraclim)
    
    # Merge data
    data = chirps.merge(era5, on = ['id','period','year','month'], how = 'outer')
    data = data.merge(terraclim, on = ['id','period','year','month'], how = 'outer')
    
    # Add column with the difference between precipitation and PET
    data['prec-pet'] = data['mean_prec'] - data['mean_pet']
    
    
    # ------------------- Calculating yearly averages ------------------------
    

    
    # First, we calculate the easy summaries
    sum1_yearly = data.drop(columns = 'month').groupby(by = ['id','period','year']).agg(
        {'mean_prec'      : 'sum',
         'mean_pet'       : 'mean',
         'mean_aet'       : 'mean',
         'mean_temp'      : 'mean', 
         'mean_rh'        : 'mean', 
         'mean_vpd'       : 'mean', 
         'prec-pet'       : lambda x: int(len(x[x<0]))
         }).reset_index()
    # Rename the length of the dry season
    sum1_yearly = sum1_yearly.rename(columns = {'prec-pet': 'ds_len'})
    
    # Then, summarise to period
    sum1_period = sum1_yearly.drop(columns = ['year']).groupby(by = ['id', 'period']).mean().reset_index()
    # print(sum1_period)
    
    # Then, calculate the dry season summaries - we just do the same, 
    # but only for those months in year when prec < pet
    # First, we calculate the easy summaries
    sum2_yearly = data[data['prec-pet'] < 0].drop(columns = ['month', 'prec-pet']).groupby(by = ['id','period','year']).agg(
        {'mean_prec'      : 'sum',
         'mean_pet'       : 'mean',
         'mean_aet'       : 'mean',
         'mean_temp'      : 'mean', 
         'mean_rh'        : 'mean', 
         'mean_vpd'       : 'mean'
         }).reset_index()
    # Rename the length of the dry season
    sum2_yearly = sum2_yearly.rename(columns = {'mean_prec': 'mean_prec_ds',
                                                'mean_pet' : 'mean_pet_ds',
                                                'mean_aet' : 'mean_aet_ds',
                                                'mean_temp': 'mean_temp_ds',
                                                'mean_rh'  : 'mean_rh_ds',
                                                'mean_vpd' : 'mean_vpd_ds',})
    
    # Then, summarise the dry-season variables to period
    sum2_period = sum2_yearly.drop(columns = ['year']).groupby(by = ['id', 'period']).mean().reset_index()
    # print(sum2_period)
    
    
    # Merging total and dry season summaries
    sum_yearly = sum1_yearly.merge(sum2_yearly, on = ['id','period','year'], how = 'outer')
    sum_period = sum1_period.merge(sum2_period, on = ['id','period'], how = 'outer')
    

    
    print('Working on CWD calculations per cell...')
    # Creating empty column to store the MCWD per period
    sum_period['mean_MCWD_aet'] = 0
    # Creating empty column to store the month peak of the dry season 
    # (month corresponding to minimum CWD)
    sum_period['mean_peak_ds'] = 0
    
    for cell in tqdm(data['id'].unique()):
        
        # Working on a period at a time
        for p in data['period'].unique():
        
            # Subsetting cell data
            cell_period_data = data[(data['id'] == cell) & (data['period'] == p)]
            cell_period_data = cell_period_data[['id','period','year','month','mean_prec','mean_aet','prec-pet']]
        
            # Making sure we have the data sorted by data
            cell_period_data = cell_period_data.sort_values(by = ['year', 'month']).reset_index(drop = True)
        
            # Getting wettest month for the cell
            wettest_month_cell = identify_wettestmonth(cell_period_data)

            
            # Calculating CWD for the cell's time series using Potential ET
            cell_period_cwd = CWD_timeseries(cell_period_data, wettest_month_cell, 'aet')
            # print(cell_period_cwd[['id','period','year','month','prec-pet', 'CWD_aet']].tail(50))
            # print(cell_period_cwd.columns)
            
            # Getting the smallest CWD per year and the corresponding month
            cell_period_mcwd = MCWD_timeseries(cell_period_cwd, wettest_month_cell)
            # print(cell_period_mcwd.head(40))
            
            # The peak of the dry season is the month with more extreme CWD
            peak_ds = cell_period_mcwd['month'].mean()
            
            # # Formatting and calculating the average MCWD in the period
            cell_period_mcwd = cell_period_mcwd['MCWD_aet'].mean()
            
            # Storing the values
            sum_period.loc[(sum_period['id'] == cell) & (sum_period['period'] == p), 'mean_MCWD_aet'] = cell_period_mcwd
            sum_period.loc[(sum_period['id'] == cell) & (sum_period['period'] == p), 'mean_peak_ds'] = peak_ds
            

    # print(sum_period)
    
    # Saving the data
    sum_yearly.to_csv(output_fp_yearly, index = False)
    sum_period.to_csv(output_fp_period, index = False)
    # 1. Merge all data together with proper suffixes
    # 2. Calculate prec chirps - pet terraclimate
    # 3. Summarise temp, prec, RH, vpd, pet, aet to year and to period
    # 4. Calculate for each year the average/total value of temp, prec, RH, vpd, pet, aet in dry season (perc < pet), then take the average for the period
    # 5. Calculate for each period the average precipitation and the average pet, idnetify the wettest month, calculate the CWD and the MCWD average in the period
    

    

    
if __name__ == '__main__':
    main()

