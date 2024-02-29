# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 17:34:42 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that gets the population by municipality data for the years 1992 to 1995
and prepares it for population density correspondence. The issue is that there 
may have been a change of municipalities between 1995 and 1996 in such a way that 
some municipalities in the period pre-1996 were split (or some were just created).
Hence, not all municipalities that appear in the shapefile of the year 2000 have
population data in the files 1992-1995.

Probably, this means that the population value for the municipalities that were
split should be divided among the new municipalities. however, since we do not 
have the look-up table with the correspondence between the old codes and the new
ones, we have to make some assumptions and approximations.

What we do is we assume the list of municiaplities in the shapefile 2000 to be
the list of municipalities for the period 1992-1995. We fill any gaps (municipalities
with no data in this early period), with the average of the municipalities 
surrounding it.

"""

"""
Imports
"""
import pandas as pd
import geopandas as gpd

import sys
sys.path.append('../utils/')
from municipalities_years_to_shp import municipality_csv_to_shp
corr = municipality_csv_to_shp()
"""
Functions
"""



"""
MAIN
"""
def main():
    
    # -------------------------------------------------------------------------
    # User inputs
    
    # The user must declare what years it want to work on
    # Start and end years
    start_year, end_year = 1985, 2012
    
    # Declare input files
    csv_fp = '../../data/processed/BRA_population/municipalities_pop_{}.csv'
    shp_fp = '../../data/shapes/municipalities_shp/municipalities_{}.shp'
    
    # Output files
    out_fp = '../../data/processed/BRA_population/municipalities_pop_{}.csv'
    
    # -------------------------------------------------------------------------
    # Calculating population density per municipality and year
    
    
    
    
    # Working on a year at a time
    for year in range(start_year, end_year + 1):
        
        print('')
        print('---------------------------------')
        print('Working on year {}'.format(year))
        
        # Read the shapefile
        shp = gpd.read_file(shp_fp.format(corr[year]))
        # Select relevant columns
        shp = shp[['state', 's_code', 'm_code', 'm_name', 'geometry']]
        
        # Creat year column in shapefile with current year
        shp['year'] = year
        
        # Read data for year
        data_y = pd.read_csv(csv_fp.format(year))
        
        # Select relevant columns only
        data_y = data_y[['year', 'm_code', 'population', 'area_km2','pop_density']]
        
        # Add the geographic information
        data_y = shp.merge(data_y, on = ['year','m_code'], how = 'left')
        
        # Municipalities in 2000 with missing data
        munic_missing = data_y[data_y['pop_density'].isna()]
        
        # Fill the empty municipalities with the average population value form the 
        # geographic neighbourhood
        while len(munic_missing) > 0:
            
            # For each missing data point, fill with average of neighbours
            for m in munic_missing['m_code']:
                
                # Gemoetry of the missing municipality
                m_geom = munic_missing[munic_missing['m_code'] == m]['geometry'].item()
                
                # Find mean population of neighbours
                neigh_pop = data_y[data_y['geometry'].touches(m_geom) & ~data_y['pop_density'].isna()]['pop_density'].mean()
                
                # Assign value to municipalities' population density
                data_y.loc[data_y['m_code'] == m, 'pop_density'] = neigh_pop
                
            # Update the list of municipalities with missing data (there may be 
            # municipalities that were empty and all the neighbours were empty too)
            munic_missing = data_y[data_y['pop_density'].isna()]
            
        # Selecting and ordering columns
        data_y = data_y[['year', 'state', 's_code', 'm_code', 'm_name', 'population', 'area_km2','pop_density']]
        
        # Data has been processed, let's save it 
        data_y.to_csv(out_fp.format(year), index = False)
            
        
        
    
    
    
if __name__ == '__main__':
    main()