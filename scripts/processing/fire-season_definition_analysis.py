# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 14:55:59 2022

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk
"""


"""
Imports
"""
import pandas as pd
import numpy as np

import burnedArea_concentrationIndex as baCI

"""
Functions
"""
def interquantile_season(ba_date, ba_col, iq):
    # Function that calculates start, end and length of fire season, defined as
    # the period of the fire year containing the 80% interquantile of the burned area
    # ba_col is the name of the column containing the burned area data
    
    # First, create column with the cumulative percentage
    # ba_date['ba_cumPerc'] = 100 * (ba_date[ba_col].cumsum() / ba_date[ba_col].sum())
    cumPerc = 100 * (ba_date[ba_col].cumsum() / ba_date[ba_col].sum())
    # print(cumPerc)
    # Fire season start
    # First, index of date with cumulative burned area percentage closest to 10%
    # HINT: remember that pd.DataFrame.idxmin() Returns index of first occurrence of minimum over requested axis.
    start_index = cumPerc.sub((100-iq)/2).abs().idxmin()
    # Hence, date 
    start_date = ba_date.loc[start_index, 'date']
    
    
    # Fire season end
    # Then, index of date with cumulative burned area percentage closest to 90%
    end_index = cumPerc.sub(100 - (100-iq)/2).abs().idxmin()
    # Hence, date 
    end_date = ba_date.loc[end_index, 'date']
    
    # Fire season length is just the difference between the indices,
    # corrected to count both the ending and starting months
    length = end_index - start_index + 1
    
    # Returning values 
    return start_date, end_date, length


def month_threshold(ba_date, ba_col):
    # Function that calculates start, end and length of fire season, defined as
    # the period of the fire year comprising all those months with burned areas
    # >= 10% of the total yearly burned area
    
    # First, calculate percentage burned area per month
    perc = 100 * ba_date[ba_col] / ba_date[ba_col].sum()
    
    # Then, subset dates with perc >= 10 %
    # Get corresponding subset of indices
    # print(perc)
    perc = perc[perc >= 10.0].reset_index()
    # print(perc)
    
    # Startint and ending months of the fire season correspond to the first 
    # and last indices
    start_index = perc.loc[0, 'index']
    start_date = ba_date.loc[start_index, 'date']
    
    end_index = perc.loc[len(perc) - 1, 'index']
    end_date = ba_date.loc[end_index, 'date']
    
    # Fire season length is just the difference between the indices,
    # corrected to count both the ending and starting months
    length = end_index - start_index + 1
    
    # Returning values 
    return start_date, end_date, length


def month_threshold_average(ba_date, ba_col, ba_average):
    # Function that calculates start, end and length of fire season, defined as
    # the period of the fire year comprising all those months with burned areas
    # >= 10% of the total yearly burned area
    # print(ba_average)
    # First, calculate percentage burned area per month using the average yearly burned area
    perc = 100 * ba_date[ba_col] / ba_average
    
    # Then, subset dates with perc >= 10 %
    # Get corresponding subset of indices
    # print(perc)
    perc = perc[perc >= 10.0].reset_index()
    # print(perc)
    
    # It may be that no months have more than 10% of the average burned area,
    # in that case, we record NaN
    if len(perc) > 0:
    
        # Startint and ending months of the fire season correspond to the first 
        # and last indices
        start_index = perc.loc[0, 'index']
        start_date = ba_date.loc[start_index, 'date']
        
        end_index = perc.loc[len(perc) - 1, 'index']
        end_date = ba_date.loc[end_index, 'date']
        
        # Fire season length is just the difference between the indices,
        # corrected to count both the ending and starting months
        length = end_index - start_index + 1
        
    else:
        start_date = np.nan
        end_date = np.nan
        length = np.nan
    
    # Returning values 
    return start_date, end_date, length

def ranking(ba_date, ba_col):
    # Function that calculates the length of the fire season by ranking months in the 
    # year and counting how many months are needed to add up to 80% of total
    # burned area in the year.

    # First, sort dataframe by burned area values from highest to lowest
    ba_date_sorted = ba_date.sort_values(by = ba_col, ascending = False).reset_index(drop = True)
    # print(ba_date_sorted)
    # Then, calculate inverse cumulative percentage
    cumPerc = 100 * (ba_date_sorted[ba_col].cumsum() / ba_date_sorted[ba_col].sum())
    # print(cumPerc)
    # Get season length as the length of the subset That first adds up to 80%
    # We choose as the cutoff month, the month closest to 80% cumulatively.
    the_index = cumPerc.sub(80.0).abs().idxmin()
    season_length = the_index + 1
    # print(the_index)
    # If the difference between the percentage burned area for this month and 80%
    # is larger than 5%, then we choose the first month with percentage smaller than 80%
    if cumPerc[the_index] - 80 > 2:
        the_index = the_index - 1
        
        if the_index < 0:
            season_length = 1
            
        elif the_index >= 0:
            season_length = the_index + 1
    # If there is no such month, then we choose the first month as the only month in the
    # burn season, hence length = 1
    # print(season_length)
    
    return season_length

def ba_concentrationIndex(ba_date, ba_col):
    # Function that calculate a burned area concentration index as an adaptation 
    # of the one used in Seasonality of Precipitation in the US and mentioned in Defining pyromes and ...
    
    # Burned Area Concentration INdex and peak of season according to this method
    ba_ci, peak_ci = baCI.burnedArea_concentrationIndex(ba_date, ba_col)
    
    # peak_ci is simply the month, get date instead
    peak_ci = ba_date.loc[ba_date['month'] == peak_ci, 'date'].item()
    # print(ba_ci)
    # print(peak_ci)
    
    return ba_ci, peak_ci
    

"""
MAIN
"""
def main():
    
    
    # ---------------------------- User inputs --------------------------------
    # Grid size
    grid_size = 100
    # Units of grid_size
    grid_units = 'km'
    
    grid_id = str(grid_size) + grid_units
    
    # Time series to work on
    start_year = 1985
    end_year = 2020
    period = str(start_year) + '-' + str(end_year)
    
    # UPDATE 15/2/2022
    # File with burned area per month-year and cell
    
    # UPDATE 15/2/2022
    # Updated to work with precipitation as well
    ba_fp = '../data/processed/grid_{}/MBfogo_burnedArea_grid{}_cerrado_{}.csv'
    # ba_fp = '../data/processed/grid_{}/summary_tables/summary_era5land_{}_{}.csv'
    
    var_name = 'ba_km2' # 'mean_prec' # 'ba_km2  ba_km2_ma50
    year_name = 'fire_year' # 'prec_year' # 'fire_year'
    
    # Lookup table to calendar year-month to fire year-month
    lkup_fp = '../data/processed/grid_{}/lookup_tables/lkup_fire-year_timeseries_cell_{}.csv'.format(grid_id, period)
    
    # File to save generated fire-season definitions data
    output_fp = '../data/processed/grid_{}/fire-season_definition/fire-season_smoothedData50_grid{}_cell_{}_2iq.csv'.format(grid_id, grid_id, period)
    # output_fp = '../data/processed/grid_{}/fire-season_definition/prec-season_data_grid{}_cell_{}.csv'.format(grid_id, grid_id, period)
    
    
    # ----------------------- Reading input -----------------------------------
    
    # ba = pd.read_csv(ba_fp)
    # Read data
    ba = []
    for year in range(start_year, end_year + 1):
        # print(year)
        ba_year = pd.read_csv(ba_fp.format(grid_id, grid_id, year))
        
        if 'year' not in ba_year.columns:
            ba_year['year'] = year
        
        ba.append(ba_year)
    ba = pd.concat(ba, ignore_index = True)
    
    # ba = pd.read_csv('../data/processed/grid_100km/smooth_timeSeries/MBfogo_burnedArea_grid100km_cerrado_smoothed_1985-2020.csv')
    
    # Reshaping from long to wide
    ba['id'] = ba['id'].apply(int).apply(str)
    ba = ba.pivot(index = ['year', 'month'], columns = 'id', values = var_name).reset_index()
    print(ba.head())
    
    lkup = pd.read_csv(lkup_fp)
    # Discarding rows with NaNs
    lkup = lkup[~lkup[year_name].isnull()].reset_index(drop = True)
    
    # Creating a column with calendar date in both dataframes
    ba['day'] = 1
    ba['date'] = pd.to_datetime(ba[['year', 'month','day']])
    lkup['day'] = 1
    lkup['date'] = pd.to_datetime(lkup[['year', 'month', 'day']])
    # erasing day columns
    ba = ba.drop(columns = 'day')
    lkup = lkup.drop(columns = 'day')

    # TODO: for third definition calcualte here the average yearly burned area per cell here to be used later
    # Calculating average yeary burned area over the time series for each cell
    # TODO; move this calculation to a function
    # First, cast ba data from wide to long format
    average_ba = ba.melt(id_vars = ['year', 'month', 'date'], var_name = 'id', value_name = var_name)
    average_ba['id'] = average_ba['id'].apply(int)
    # Merging to lookup table
    average_ba = average_ba.merge(lkup, on = ['year', 'month', 'date', 'id'], how = 'right')
    
    # Now we can first calcualte total burned area per fire year
    average_ba = average_ba[['id',year_name, var_name]].groupby(by = ['id', year_name]).sum().reset_index()
    # And now we calculate average yearly burned area per grid cell
    average_ba = average_ba[['id', var_name]].groupby(by = ['id']).mean().reset_index()
    # print(average_ba)

    # ---------------- Fire season calculation --------------------------------
    
    # First, extracting list of fire year labels 
    fire_years = lkup[year_name].unique()
    
    # List of cells
    cells_list = lkup['id'].unique()
    # print(cells_list[0:3])
    
    # Create list to accumulate dataframes (one per grid cell and fire year combination)
    list_df = []
    
    # Working on each cell at a time 
    # TODO: write possibility to do this in parallel
    for cell in cells_list:
        
        cell = 440
    
        print('\n Working on cell {}'.format(cell))    
    
        # Calculate seasonality per fire year
        for fyear in fire_years:
    
            # print('Working on fire year {}'.format(fyear))
            # Subset calendar dates belonging to this fire year
            lkup_fyear = lkup.loc[(lkup[year_name] == fyear) & (lkup['id'] == cell), 'date']
            # print(lkup_fyear)
            print(lkup_fyear)
            
            # print(ba.head())
            # print(ba.columns)
            # print(ba[12])
            
            # print(ba.loc[ba['date'].isin(lkup_fyear), ['month', 'date', '12']])
            
            # Use these dates to subset burned area data belonging to this fire year and cell
            ba_fyear = ba.loc[ba['date'].isin(lkup_fyear), ['month', 'date', str(int(cell))]].reset_index(drop = True)
            print(ba_fyear)
        #     break
        # break
            # It might be that there are no fires in the fire year, then we fill the year with NaNs
            if ba_fyear[str(int(cell))].sum() != 0.0:
                # Obtaining peak of fire season
                peak = ba_fyear.loc[ba_fyear[str(int(cell))].idxmax(), 'date']
                
                # Create data frame to store seasonality information for this fire year and cell
                cell_fyear_season = pd.DataFrame({'id': cell, 'fyear': fyear, 'peak_season': [peak]})
                
                # Fire sesaon per definition --------------------------------------
                
                # DEFINITION 1: using an 80% interquantile
                cell_fyear_season['iq80_st'], cell_fyear_season['iq80_end'], cell_fyear_season['iq80_length'] = interquantile_season(ba_fyear, str(int(cell)), 80)
                
                # DEFINITION 2: using a monthly thresold over the yearly burned area
                cell_fyear_season['monthThresh_st'], cell_fyear_season['monthThresh_end'], cell_fyear_season['monthThresh_length'] = month_threshold(ba_fyear, str(int(cell)))
                
                # DEFINITION 3: using a monthly thershold over the average year burned area
                cell_fyear_season['monthThreshAver_st'], cell_fyear_season['monthThreshAver_end'], cell_fyear_season['monthThreshAver_length'] = month_threshold_average(ba_fyear, str(int(cell)), average_ba.loc[average_ba['id'] == cell, var_name].item())
                
                # DEFINITION 4: length of season by ranking months according to burned area
                cell_fyear_season['ranking_length'] = ranking(ba_fyear, str(int(cell)))
                
                # DEFINITION 5: Adaptation of the rainfall concentration index
                cell_fyear_season['ba_CI'], cell_fyear_season['peak_CI'] = ba_concentrationIndex(ba_fyear, str(int(cell)))
                
                # DEFINITION 6: using an 90% interquantile
                cell_fyear_season['iq90_st'], cell_fyear_season['iq90_end'], cell_fyear_season['iq90_length'] = interquantile_season(ba_fyear, str(int(cell)), 90)
                
                # print(cell_fyear_season)
                list_df.append(cell_fyear_season)
                
            else:
                cell_fyear_season = pd.DataFrame({'id': cell, 'fyear': fyear, 'peak_season': [np.nan],
                                                  'iq_st':np.nan, 'iq_end': np.nan, 'iq_length': np.nan,
                                                  'monthThresh_st': np.nan, 'monthThresh_end': np.nan, 'monthThresh_length': np.nan,
                                                  'monthThreshAver_st':np.nan, 'monthThreshAver_end':np.nan, 'monthThreshAver_length': np.nan,
                                                  'ranking_length': np.nan,
                                                  'ba_CI': np.nan, 'peak_CI': np.nan})
                list_df.append(cell_fyear_season)
            break
        

        break
    cell_fyear_season = pd.concat(list_df, ignore_index = True)
    
    cell_fyear_season.to_csv(output_fp, index = False)
if __name__ == '__main__':
    main()

