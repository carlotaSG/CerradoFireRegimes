# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 10:21:32 2022

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that joins spatially, into a single MultiPolygon, those MBfogo polygons 
that belong to the same month and that are spatially close to one another.

These polygons are considered to be a single MultiPolygon based on 
a certain buffering distance. To test the effects of different buffers, the script
generates different shapefiles considering buffers equal to:
    - 1 to 10 Landsat pixels (30 to 300 m)
    - 500, 750 and 1,000 m.
    


"""

"""
Imports
"""
import geopandas as gpd
import pandas as pd
import dask_geopandas as dgpd
from libpysal.weights import fuzzy_contiguity
import sys
# from shapely.geometry import MultiPolygon

sys.path.append('../utils/')

from multiprocessing import Pool

from create_yearMonthID import create_yearMonthID

from pytictoc import TicToc
t = TicToc()


"""
Functions
"""

def id_connectedComponents(gdf, buff, month, last_W):
    # Function that gets a GeoDataFrame and, on a month by month basis, 
    # identifies those polygons that interset one another with a certain buffer
    
    # Creating list to accumulate GeoDataFrames with identified connected components 
    # pols_connected_id = []
    
    # If working on the first month, then starting
    if month == 1:
        # Initialising counter to keep track of the right identifier of the connected
        # component
        last_W = -1
    
    # # Operating on a monthly basis
    # for month in range(1, 12 + 1):
    t.tic()
    print('Identifying connected components for month {}'.format(month))
    
        # # Subsetting polygons mapped on month
        # pol_month = gdf[gdf['month'] == month]
    
    # Identifying the polygons that are within a distance from one another
    W = fuzzy_contiguity(gdf, buffering = True, buffer = buff, silence_warnings = True)
        
    # Correcting to accumulate index through months
    conn_ids = W.component_labels + last_W + 1
    # print(conn_ids)
        
    gdf['conn_ids'] = conn_ids
        
        # pols_connected_id.append(pol_month)
        
    # Updating the last weight position
    last_W = max(conn_ids)
        
    del W, conn_ids
        
    # # Identifying the polygons that are within a distance from one another
    # W = fuzzy_contiguity(gdf, buffering = True, buffer = buff, silence_warnings = True) 
    
    # gdf['conn_ids'] = W.component_labels
    
        
    t.toc('Connected components for month {} identified in'.format(month))
        
    # pols_connected_id = gpd.GeoDataFrame(pd.concat(pols_connected_id))
    
    return gdf, last_W


def create_lookupID(gdf):
    # Fnction that creates a lookup table that relates the original polygon IDs
    # to the new polygon IDs (the connected polygons)
    
    # Create lookup base
    lookup = gdf[['ID','year','month','conn_ids']]
    # print(lookup)
    
    # First, create a unique ID for each new ID
    # Select the new ids
    new_ids = lookup.drop(columns = ['ID']).drop_duplicates(subset = 'conn_ids').reset_index(drop = True)
    # print(new_ids)
    # Create an ID with the structure <year><month>_<dataframe_index>, where
    # dataframe_index is already the connected_id (as from previous operation)
    new_ids = create_yearMonthID(new_ids).reset_index()
    # print(new_ids)
    # Renaming connected ID column
    new_ids = new_ids.rename(columns = {'ID': 'ID_sjoin'})
    
    # Merge this new datafrme back to lookup
    lookup = lookup.merge(new_ids[['ID_sjoin','conn_ids']], on = ['conn_ids'], how = 'left')
    # print(lookup)
    
    # Freeing memory
    del new_ids
    
    # Drop the conn_ids column as no longer needed
    # lookup = lookup.drop(columns = 'conn_ids')
    
    return lookup


def dissolve_polygons(gdf):
    
    gdf = gdf.dissolve(by = 'ID_sjoin').reset_index()
    
    return gdf


def aggregate_connectedPolygons(gdf):
    # Function that aggregates same-group polygons into MultiPolygon. First, it
    # merges the connected-polygon ID using the lkup table
    
    # print(gdf)
    # print(lkup)
  
    # Merge lookup information back to polygons GeoDataFrame
    # gdf = gdf.merge(lkup, on = ['conn_ids','year','month'], how = 'left')
    # print(gdf)
    
    # and aggregate by ID_sjoin
    # TODO: to make this faster, I could separate the repeated ones (do not need dissolve)
    # and the non-repeated ones (require dissolve)
    # I could also parallelise over months
    # gdf = gdf[['ID_sjoin','year','month','geometry']].dissolve(by = 'ID_sjoin').reset_index()
    
    
    
    
    # Try using Dask
    
    # https://dask-geopandas.readthedocs.io/en/latest/guide/dissolve.html
    # accessed 25/07/2022
    
    # Perhaps I could still do better and take out those polygons that are not part
    # of a connected component
    # First, getting list of repeated rows
    # mask = gdf.ID_sjoin.duplicated(keep = False)
    
    # Hence, non-repeated elements
    # gdf_nonconnected = gdf.loc[~mask, ['ID_sjoin','year','month','geometry']]
    
    # And repeated elements
    # gdf_connected = gdf.loc[mask, ['ID_sjoin','year','month','geometry']]
    
    # del gdf
    
    # gdf = gdf[['ID_sjoin','year','month','geometry']]
    
    ddf = dgpd.from_geopandas(gdf, npartitions=35)
    
    shuffled = dissolve_shuffle(ddf, 'conn_ids')
    
    gdf = shuffled.compute()
    
    # gdf = gpd.GeoDataFrame(pd.concat([gdf_nonconnected, gdf_connected], ignore_index = True), crs = 'epsg:5880')
    
    # del gdf_connected, gdf_nonconnected
    
    # # First, getting list of repeated rows
    # mask = gdf.ID_sjoin.duplicated(keep = False)
    
    # # Hence, non-repeated elements
    # gdf_nonConnected = gdf.loc[~mask, ['ID_sjoin','year','month','geometry']]
    
    # # And repeated elements
    # gdf_connected = gdf.loc[mask, ['ID_sjoin','year','month','geometry']]
    
    
    # # Freeing memory
    # del gdf
    
    # # Now, for those repeated elements, we are going to parallelise over months,
    # # so that each worker deals with a different month, in parallel
    # # This takes advantage of the fact that the connected parts belong to a unique month
    # poolInput = [gdf_connected[gdf_connected['month'] == month] for month in gdf_connected['month'].unique()]
    
    # # print(poolInput)
    
    # # Freeing memory
    # del gdf_connected
    
    # # Parallelising dissolve function
    # with Pool(8) as p:
    #     poolOutput = p.map(dissolve_polygons, poolInput)
        
    # # Output is a list of GeoDataframes - Now we have the dissolved polygons, of those that were connected
    # poolOutput = gpd.GeoDataFrame(pd.concat(poolOutput, ignore_index = True))
    
    # # Now we concatenate this list of MultiPolygons to those polygons that were not
    # # part of a group
    # gdf = gpd.GeoDataFrame(pd.concat([gdf_nonConnected, poolOutput], ignore_index = True))
    
    # # Freeing memory
    # del poolOutput, gdf_nonConnected
    
    
        
    
    
        
    
    # # Using dissolve to create MultiPolygons from Polygons may be creating issues due to 
    # # floating point precision inaccuracies. Hence, all of a sudden polygons joined
    # # are not maintaining the shape they should be keeping. This is because dissolve
    # # uses a unary_union method, instead of just plugging the polygons into a 
    # # list of polygons (hence creating a MultiPolygon in a simple manner).
    
    # # Therefore, I create the Multipolygons "manually"
    
    # # Select relevant set of columns
    # gdf = gdf[['ID_sjoin','year','month','geometry']]
    
    # # First, create GeoDataframe with the unique connected id information, the 
    # # year and month, with a dummy geometry column (non-correct geometries)
    # gdf_connected = gdf.drop_duplicates(subset = ['ID_sjoin','year','month']).reset_index(drop = True)    
    
    # # print(len(gdf))
    # # print(len(gdf_connected))
    
    # for con_id in gdf_connected['ID_sjoin']:
        
    #     # If there is more than one polygon belonging to this group, then we 
    #     # cast the various geometries into a MultiPolygon
    #     if len(gdf[gdf['ID_sjoin'] == con_id]) > 1:
            
    #         # print('here')
            
    #         list_polygons = gdf.loc[gdf['ID_sjoin'] == con_id, 'geometry'].to_list()
    #         # print(list_polygons)
            
    #         list_polygons = to_polygons(list_polygons)
    #         # print(list_polygons)
            
    #         multipol = MultiPolygon(list_polygons)
    #         # print(multipol)
            
    #         # print(multipol.wkt)
            
    #         # print(gpd.GeoSeries([multipol], crs = 'epsg:5880'))
            
    #         # print(multipol)
    #         # print(multipol)
    #         # print(gdf_connected.loc[gdf_connected['ID_sjoin'] == con_id, :])
            
    #         # print(gpd.GeoSeries([multipol]))
    #         # print(len(gpd.GeoSeries([multipol])))
            
    #         gdf_connected.loc[gdf_connected['ID_sjoin'] == con_id, 'geometry'] = gpd.GeoSeries([multipol], crs = 'epsg:5880')
            
    #         # gdf_connected.loc[gdf_connected['ID_sjoin'] == con_id, 'geometry'] = shapely.geometry.MultiPolygon(gdf.loc[gdf['ID_sjoin'] == con_id, 'geometry'].to_list())
    
    # return gdf_connected
    
    return gdf



# def to_polygons(geometries):
#     # Generator function that makes a list of polygons out of a list of both
#     # polygons and MultiPolygons
#     # From https://stackoverflow.com/questions/60691287/convert-list-of-multipolygons-and-polygons-to-a-single-multipolygon
#     # Accessed on 25/07/2022
#     for geometry in geometries:
#         if isinstance(geometry, Polygon):
#             yield geometry
#         else:
#             yield from geometry


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
    
    # Asking the user whether to specify the year to work on:
    year = input("Input year to work on:\n")
    
    # If incorrect year input, notify and exit program
    if len(year) != 4:
        print('\nWARNING: wrong year input.\nExiting program.')
        sys.exit(1)
    else:
        year = int(year)
    
    # grid_id = '32deg'
    # cell_id = 25
    
    # Shapefile of polygons
    # pol_fp = '../data/testing_buffers_spatial/MBfogo_c10_cell_{}_id{}_polygons_annual.shp'
    pol_fp = '../../data/tests/fire_polygons/individual_noMonths/MBfogo_cerrado_polygons_{}_fixedGeom.shp'.format(year)
    
    # Output files
    # csv_out_fp = '../data/processed/MBfogo_c10_processed_/MBfogo_c10_cell_{}_id{}_polygons_2019_buff{}_lookup.csv'
    # shp_out_fp = '../data/testing_buffers_spatial/MBfogo_c10_cell_{}_id{}_polygons_annual_buff{}.shp'
    shp_out_fp = '../../data/tests/fire_polygons/buff{}/MBfogo_cerrado_polygons_{}_fixedGeom_buff{}.shp'
    
    # List of buffers to work on 
    list_buffers = [15, 45, 75, 150, 250, 500]#, 250, 500]
    
    # # For the annual version, we only merge those polygons that are touching one another
    # list_buffers = [5, 15]
    
    
    
    # Working on a buffer at a time
    for buffer in list_buffers:
        
        print('Working on buffer {}'.format(buffer))
        
        # Initialising id counter
        last_W = -1
        
        
        # Working on a month by month basis
        for month in [1]:
            
            print('Working on month {}'.format(month))
            
            # Reading and subsetting relevant data
            pol = gpd.read_file(pol_fp)
            
            # Renaming raster_val column to month, all of them will have month 1
            pol = pol.rename(columns = {'raster_val':'month'})
            # Creating column with year
            pol['year'] = year
            
            # Subsetting monthly data
            pol = pol[pol['month'] == month]

            # Converting to metre CRS
            print("Converting polygon's CRS")
            crs_original = pol.crs
            pol = pol.to_crs('epsg:5880')
        
        
            
        
            # Identifying the connected components
            t.tic()
            print('Identifying connected polygons per month...')
            connected_pol, last_W = id_connectedComponents(pol, buffer, month, last_W)
            t.toc('Connected polygons identified in')
            
            # print(connected_pol.head())
            
            del pol
        
            # Here, we should separate those polygons that are to be dissolved, 
            # from the ones that are not. And from these last ones, we can already write
            # to shapefile those corresponding to month-1
            print('Sorting non-connected components...')
            t.tic()
            # 1. Finding the duplicates
            duplicates = connected_pol['conn_ids'].duplicated(keep = False)
            # 2. Hence, separating those polygons not to be dissolved, and those that are
            pols_notDissolve = connected_pol[~duplicates]
            connected_pol = connected_pol[duplicates]
            
            
            print('Saving non-connected polygons for month {}...'.format(month))
            # Further, save and remove those polygons not to be dissolved from month-1
            if month == 1:
                # If we are working on the first round of month pairs, create output file and save
                pols_notDissolve[['conn_ids','year','month','geometry']].to_crs(crs_original).to_file(shp_out_fp.format(buffer, year, buffer))
            else:
                # Write these polygons appending them to output file
                pols_notDissolve[['conn_ids','year','month','geometry']].to_crs(crs_original).to_file(shp_out_fp.format(buffer, year, buffer), mode = 'a')
            # Freeing memory
            del pols_notDissolve
            t.toc('Non-connected sorted in')
        

            # print(len(connected_pol))
            # If there are any polygons to connect, aggregate to a MultiPolygon
            if len(connected_pol) > 0:
                t.tic()
                print('Dissolving polygons into MultiPolygons...')
                connected_pol = aggregate_connectedPolygons(connected_pol)
                t.toc('Polygons dissolved in')
            
                # Now, save the shapefile with the connected polygons
                t.tic()
                print('Saving files for month {}'.format(month))
                # TODO: convert back to initial CRS!!
                connected_pol = connected_pol.to_crs(crs_original)
                connected_pol[['conn_ids','year','month','geometry']].to_file(shp_out_fp.format(buffer, year, buffer), mode = 'a')
                t.toc('Files saved in')
            
            del connected_pol
            
            print('')
            
        
        print('Finished merging polygons of the same month.')
        
        print('\n \n')

        # Creating unique ID for each polygon
        t.tic()
        print('Creating unique polygon IDs...')
        # First, read polygons again
        pols = gpd.read_file(shp_out_fp.format(buffer, year, buffer))
        # Then, order by year and month
        pols = pols.sort_values(by = ['year','month']).reset_index(drop = True)
        # Adding column with the IDs
        pols = create_yearMonthID(pols)
        
        # print(pols.head())
        # Modifying column name for clarity
        pols = pols.rename(columns = {'ID':'ID_sjoin'})
        # Saving again polygons with the right IDs
        pols.drop(columns = {'conn_ids'}).to_file(shp_out_fp.format(buffer, year, buffer))
        t.toc('IDs table created and saved in')
        print('\n \n \n ----------------------------------------------')
        

        
    # Setting on slice DataFrame warning
    pd.options.mode.chained_assignment = 'warn'

if __name__ == '__main__':
    main()

