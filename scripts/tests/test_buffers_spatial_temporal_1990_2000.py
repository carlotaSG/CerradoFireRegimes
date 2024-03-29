# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 11:55:47 2022

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that joins spatially, into a single MultiPolygon, those MBfogo polygons 
that belong to consecutive months and that are spatially close to one another.

These polygons are considered to be a single MultiPolygon based on 
a certain buffering distance. Less than one Landsat pixel in this case.
    

"""


"""
Imports
"""
import geopandas as gpd
import pandas as pd
import dask_geopandas as dgpd
from libpysal.weights import fuzzy_contiguity
# from shapely.geometry import MultiPolygon
import sys
import os

from multiprocessing import Pool

sys.path.append('../utils/')
from create_yearMonthID import create_yearMonthID

from pytictoc import TicToc
t = TicToc()


"""
Functions
"""

def id_connectedComponents(gdf, buff):
    # Function that gets a GeoDataFrame with polygons for two consecutive months,
    # identifies those polygons that interset one another with a certain buffer
    
    # Identifying the polygons that are touching one another with a 10 m buffer
    # (5 m from each polygon)
    W = fuzzy_contiguity(gdf, buffering = True, buffer = 5, silence_warnings = True)
        
        
    gdf['connected_ids'] = W.component_labels

    
    return gdf





def aggregate_connectedPolygons(gdf):
    # Function that aggregates same-group polygons into MultiPolygon. First, it
    # merges the connected-polygon ID using the lkup table

    
    ddf = dgpd.from_geopandas(gdf, npartitions=30)
    
    shuffled = dissolve_shuffle(ddf, 'connected_ids', aggfunc={
        'conn_ids' : 'first',
        'year' : 'first',
         "month" : "max",
         # "pop_est": ["min", "max"],
     })
    
    gdf = shuffled.compute()
    
    # Getting rid of any left inner lines
    gdf['geometry'] = gdf['geometry'].buffer(0.0001, join_style = 2)
    # Finally, shrink back geometries again
    gdf['geometry'] = gdf['geometry'].buffer(-0.0001, join_style = 2)    
    
    return gdf



def dissolve_shuffle(ddf, by=None, **kwargs):
    """Shuffle and map partition"""

    meta = ddf._meta.dissolve(by=by, as_index=False, **kwargs)

    shuffled = ddf.shuffle(
        by, npartitions=ddf.npartitions, shuffle="tasks", ignore_index=True
    )

    return shuffled.map_partitions(
        gpd.GeoDataFrame.dissolve, by=by, as_index=False, meta=meta, **kwargs
    )
        

"""
MAIN
"""
def main():
    
    # Setting off slice DataFrame warning
    pd.options.mode.chained_assignment = None  # default='warn
    
    # User inputs ------------------------------------------------------------
    
    # Asking the user to specify the year to work on:
    year = input("Input year to work on:\n")
    
    # If incorrect year input, notify and exit program
    if len(year) != 4:
        print('\nWARNING: wrong year input.\nExiting program.')
        sys.exit(1)
    else:
        year = int(year)
    
    # Input file
    pol_fp = '../../data/tests/fire_polygons/buff{}/MBfogo_c10_cerrado_polygons_subset_{}_buff{}.shp'
    
    
    # Output directory
    out_dir = '../../data/tests/fire_polygons/buff{}_consecutiveMonths/'
    # Output file
    shp_out_fp = out_dir + 'MBfogo_c10_cerrado_polygons_subset_{}_buff{}.shp'

    
    # List of buffers to work on 
    # they correspond to distances of: 30, 90, 150, 300, 500, 1000 metres
    list_buffers = [500]#15, 45, 75, 150, 250, 
    
    # grid_id = '32deg'
    # cell_id = 25
    
    # The month we wish to start working on
    start_month = 1
    # The last month we want to work on
    end_month = 12
    
    # Reading file -----------------------------------------------------------
    
    
    
    # Working on a buffer at a time - this is only relevant to choose the input file to work on
    for buffer in list_buffers:
        
        print('Working on buffer {}'.format(buffer))
        
        # Reading polygons for first month only
        pols = gpd.read_file(pol_fp.format(buffer, year, buffer))
        pols = pols[pols['month'] == start_month]
        pols = pols[['conn_ids','year','month','geometry']]
        

        # Converting to metre CRS
        print("Converting polygon's CRS")
        crs_original = pols.crs
        pols = pols.to_crs('epsg:5880')
        
        # Creating output directory if it doesn't exist
        if not os.path.exists(out_dir.format(buffer)):
            os.makedirs(out_dir.format(buffer))
        
        
        
        # We work on one month's data at a time. Merging it with the previous month
        # Starting on second month, since start_month hos no prior month
        # to merge it to
        for month in range(start_month + 1, end_month + 1):
            
            # Reading this month's data and appending it to the previous month's
            # polygons
            pol_month = gpd.read_file(pol_fp.format(buffer, year, buffer), include_fields=["conn_ids", "year", "month"])
            pol_month = pol_month[['conn_ids','year','month','geometry']]
            pol_month = pol_month[pol_month['month'] == month]
            pol_month = pol_month.to_crs('epsg:5880')
            pols = gpd.GeoDataFrame(pd.concat([pols, pol_month], ignore_index = True), crs = 'epsg:5880')
            
            # Freeing memory
            del pol_month
            
            # Identifying the connected components
            t.tic()
            print('Identifying connected polygons for months {} and {}...'.format(month, month - 1))
            pols = id_connectedComponents(pols, buffer)
            t.toc('Connected polygons identified in')
            
            # Here, we should separate those polygons that are to be dissolved, 
            # from the ones that are not. And from these last ones, we can already write
            # to shapefile those corresponding to month-1
            print('Sorting non-connected components...')
            t.tic()
            # 1. Finding the duplicates
            duplicates = pols['connected_ids'].duplicated(keep = False)
            # 2. Hence, separating those polygons not to be dissolved, and those that are
            pols_notDissolve = pols[~duplicates]
            pols = pols[duplicates]
            
            print('Saving non-connected polygons from previous month (month = {})'.format(month - 1))
            # Further, save and remove those polygons not to be dissolved from month-1
            if month - 1 == start_month:
                if not pols_notDissolve.loc[pols_notDissolve['month'] == month - 1, ['conn_ids','year','month','geometry']].empty:
                    # If we are working on the first round of month pairs, create output file and save
                    pols_notDissolve.loc[pols_notDissolve['month'] == month - 1, ['conn_ids','year','month','geometry']].to_crs(crs_original).to_file(shp_out_fp.format(buffer, year, buffer))
            else:
                if not pols_notDissolve.loc[pols_notDissolve['month'] == month - 1, ['conn_ids','year','month','geometry']].empty:
                    # Write these polygons appending them to output file
                    pols_notDissolve.loc[pols_notDissolve['month'] == month - 1, ['conn_ids','year','month','geometry']].to_crs(crs_original).to_file(shp_out_fp.format(buffer, year, buffer), mode = 'a')
            # Removing these polygons
            pols_notDissolve = pols_notDissolve[pols_notDissolve['month'] == month]
            t.toc('Non-connected sorted in')
            
            # 3. Dissolve connected polygons to a MultiPolygon
            t.tic()
            print('Dissolving polygons from month {} into MultiPolygons...'.format(month))
            pols = aggregate_connectedPolygons(pols)
            t.toc('Polygons dissolved in')
            
            # Finally, concatenating back these two set of polygons to see if they
            # are connected to any next month's polygons
            pols = gpd.GeoDataFrame(pd.concat([pols_notDissolve, pols], ignore_index = True), crs = 'epsg:5880')
            pols = pols[['conn_ids','year','month','geometry']]
            
            
        
        print('Saving last polygons to file')
        if not pols.empty:
            # Saving the last polygons, appending them to file
            pols.to_crs(crs_original).to_file(shp_out_fp.format(buffer, year, buffer), mode = 'a')
                
            
        
        

        
        
    
        
    # Setting on slice DataFrame warning
    pd.options.mode.chained_assignment = 'warn'

if __name__ == '__main__':
    main()
    
    


#==============================================================================
# UPDATE 22/02/2023

# Archiving here only for legacy purposes chunks of code that have
# been commented out


# User inputs ------------------------------------------------------------
# # List of buffers of the polygons months output by test_buffers_spatial
# list_buffers = [15, 45, 75]#, 150, 250, 500]

# # Shapefile of polygons
# pol_fp = '../data/testing_buffers_spatial/MBfogo_c10_cell_{}_id{}_polygons_2019_buff{}.shp'

# # Output files
# # csv_out_fp = '../data/testing_buffers_spatial/MBfogo_c10_cell_{}_id{}_polygons_2019_buff{}_lookup.csv'
# shp_out_fp = '../data/testing_buffers_spatial_temporal/MBfogo_c10_cell_{}_id{}_polygons_2019_buff{}_consecutiveMonths.shp'

        # del pol
        
        # # Creating lookup table between the original polygon IDs and the connected IDs
        # t.tic()
        # print('Creating lookup table between original and connected IDs...')
        # lookup = create_lookupID(connected_pol)
        # t.toc('Lookup table created in')

# def create_lookupID(gdf):
#     # Function that creates a lookup table that relates the original polygon IDs
#     # to the new polygon IDs (the connected polygons)
    
#     # Create lookup base
#     lookup = gdf[['ID','year','month','connected_ids']]
#     # print(lookup)
    
#     # First, create a unique ID for each new ID
#     # Select the new ids
#     new_ids = lookup.drop(columns = ['ID']).drop_duplicates(subset = 'connected_ids').reset_index(drop = True)
#     # print(new_ids)
#     # Create an ID with the structure <year><month>_<dataframe_index>, where
#     # dataframe_index is already the connected_id (as from previous operation)
#     new_ids = create_yearMonthID(new_ids).reset_index()
#     # print(new_ids)
#     # Renaming connected ID column
#     new_ids = new_ids.rename(columns = {'ID': 'conn_ids'})
    
#     # Merge this new datafrme back to lookup
#     lookup = lookup.merge(new_ids[['conn_ids','connected_ids']], on = ['connected_ids'], how = 'left')
#     # print(lookup)
    
#     # Freeing memory
#     del new_ids
    
#     # Drop the connected_ids column as no longer needed
#     # lookup = lookup.drop(columns = 'connected_ids')
    
#     return lookup


# def dissolve_polygons(gdf):
    
#     gdf = gdf.dissolve(by = 'conn_ids').reset_index()
    
#     return gdf

# def id_connectedComponents(gdf, buff):
    # Function that gets a GeoDataFrame with polygons for two consecutive months,
    # identifies those polygons that interset one another with a certain buffer
    
    # # Creating list to accumulate GeoDataFrames with identified connected components 
    # pols_connected_id = []
    # # Initialising counter to keep track of the right identifier of the connected
    # # component
    # last_W = -1
    
    # # Operating on a monthly basis
    # for month in range(1, 12 + 1):
        # t.tic()
        # print('Identifying connected components for month {}'.format(month))
    
        # # Subsetting polygons mapped on month
        # pol_month = gdf[gdf['month'] == month]
        
        
# def aggregate_connectedPolygons(gdf):
    # Function that aggregates same-group polygons into MultiPolygon. First, it
    # merges the connected-polygon ID using the lkup table
    
    # print(gdf)
    # print(lkup)
  
    # Merge lookup information back to polygons GeoDataFrame
    # gdf = gdf.merge(lkup, on = ['connected_ids','year','month'], how = 'left')
    # print(gdf)
    
    # and aggregate by conn_ids
    # TODO: to make this faster, I could separate the repeated ones (do not need dissolve)
    # and the non-repeated ones (require dissolve)
    # I could also parallelise over months
    # gdf = gdf[['conn_ids','year','month','geometry']].dissolve(by = 'conn_ids').reset_index()
    
    
    
    
    # Try using Dask
    
    # https://dask-geopandas.readthedocs.io/en/latest/guide/dissolve.html
    # accessed 25/07/2022
    
    # Perhaps I could still do better and take out those polygons that are not part
    # of a connected component
    # First, getting list of repeated rows
    # mask = gdf.conn_ids.duplicated(keep = False)
    
    # Hence, non-repeated elements
    # gdf_nonconnected = gdf.loc[~mask, ['conn_ids','year','month','geometry']]
    
    # And repeated elements
    # gdf_connected = gdf.loc[mask, ['conn_ids','year','month','geometry']]
    
    # del gdf
    
    # gdf = gdf[['conn_ids','year','month','geometry']]
    
# def id_connectedComponents(gdf, buff):
        
        # Correcting to accumulate index through months
        # connected_ids = W.component_labels + last_W + 1
        # print(connected_ids)
        
            
        # pols_connected_id.append(pol_month)
        
        # Updating the last weight position
        # last_W = max(connected_ids)
        
        # del pol_month, W, connected_ids
        
    # # Identifying the polygons that are within a distance from one another
    # W = fuzzy_contiguity(gdf, buffering = True, buffer = buff, silence_warnings = True) 
    
    # gdf['connected_ids'] = W.component_labels
    
        
        # t.toc('Connected components for month {} identified in'.format(month))
        
    # pols_connected_id = gpd.GeoDataFrame(pd.concat(pols_connected_id))