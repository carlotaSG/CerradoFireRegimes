# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 11:55:47 2022

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that joins into a single MultiPolygon those MBfogo polygons 
that belong to consecutive months and that are spatially touching one another
(within 10 m).

Input files:
    ../../data/temp/<year>/MBfogo_c20_cerrado_polygons_<year>_buff<buffer_label>.shp

Output files:
    ../../data/temp/<year>_CM/MBfogo_c20_cerrado_polygons_<year>.shp
    
WARNING: the input files are delete at the end of this script.

"""


"""
Imports
"""
import geopandas as gpd
import pandas as pd
import dask_geopandas as dgpd
# from libpysal.weights import fuzzy_contiguity
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix
import numpy as np
import sys
import os, shutil


sys.path.append('../utils/')

import fileList

from pytictoc import TicToc
t = TicToc()


"""
Functions
"""
def read_monthPolygons_convertCRS(gdf_fp, month_indicator, month):
    """
    Function that reads the polygons corresponding to month. It does so by reading
    only a slice of rows indicated by month_indicator (which contains the relation index-month) 
    and month. It also converts the polygons to the EPSG:5880 CRS, but it records
    the original CRS and returns it.

    Parameters
    ----------
    gdf_fp : string
        Filepath to GeoDataFrame containing the polygons to be read.
    month_indicator : DataFrame
        DataFrame with column (index) | year | month, to be used to learn the 
        slice of indices with month=month.
    month : integer
        The month of interest for which to read the polyons.

    Returns
    -------
    polygons : GeoDataFrame
        GeoDataFrame containing only those polygons from gdf_fp belonging to month.
    crs_original : string
        The original CRS of gdf_fp.

    """
    # print(month_indicator)
    
    # Get the row indices that contain start_month's data
    month_index = month_indicator[month_indicator['month'] == month].index
    
    # print(month_index)
    
    # Reading and subsetting start_month's data. 
    # As this is the first month, there is no previous month to merge it to
    polygons = gpd.read_file(gdf_fp, rows = slice(month_index[0], month_index[-1]),
                             engine = 'pyogrio')
    polygons = polygons[['conn_ids','year','month','geometry']]
    
    # Freeing memory 
    del month_index

    
    # Recording original CRS and converting to metre CRS
    print("Converting polygon's CRS")
    crs_original = polygons.crs
    polygons = polygons.to_crs('epsg:5880')
    
    return polygons, crs_original



def id_connectedComponents(gdf_prev_m, gdf_next_m, num_process):
    
    # First, reset the indices of the two dataframes
    gdf_prev_m = gdf_prev_m.reset_index(drop = True)
    gdf_next_m = gdf_next_m.reset_index(drop = True)

    
    # Creating new column identifier for each month, which has continuity from 
    # gdf_prev_m tp gdf_next_m. This column will be the index of the adjacency
    # matrix of the polygons connectivity
    gdf_prev_m['ID_prev'] = gdf_prev_m.index
    gdf_next_m['ID_next'] = gdf_next_m.index + gdf_prev_m['ID_prev'].max() + 1
    
    # Recording the last index in matrix
    last_index = gdf_next_m['ID_next'].max()
    
    # Recording the original crs
    crs_m = gdf_prev_m.crs
    
    # Creating dask-geopandas objects out of each
    dgdf_prev_m = dgpd.from_geopandas(
        gdf_prev_m,
        npartitions = num_process)
    dgdf_next_m = dgpd.from_geopandas(
        gdf_next_m,
        npartitions = num_process)
    
    # Add buffer of 5 m around each polygon using the dask infrastructure
    # so that we make sure we are considering as together all plygons at 10 m.
    # This is a sanity measure to avoid point float precission as much as possible
    dgdf_prev_m['geometry'] = dgdf_prev_m['geometry'].buffer(5, join_style = 2).compute()
    dgdf_next_m['geometry'] = dgdf_next_m['geometry'].buffer(5, join_style = 2).compute()
    
    
    # Drawing the spatial partitions for each one
    dgdf_prev_m = dgdf_prev_m.spatial_shuffle()
    dgdf_next_m = dgdf_next_m.spatial_shuffle()
    
    
    # Doing the spatial join
    # TODO: we are doing intersect if adding the buffer
    dgdf = dgpd.sjoin(dgdf_prev_m, dgdf_next_m, predicate = 'intersects')
    # dgdf = dgdf_prev_m.sjoin(dgdf_next_m, predicate = 'intersects')    
    dgdf = dgdf.compute()
    dgdf = dgdf.reset_index(drop = True)
    
    # Freeing memory
    del dgdf_prev_m, dgdf_next_m
    
    # Casting information collected into a sparse matrix, using the created IDs
    # as matrix index. We do this using scipy sparse matrix
    sparse_matrix = csr_matrix((np.ones(len(dgdf)),                            # Dummy sparse matrix values
                                (dgdf['ID_prev'], dgdf['ID_next'])),           # (row indices, column indices)
                               shape = (last_index + 1, last_index + 1))       # Shape of matrix is the total
                                                                               # number of polygons between the two months
    # Freeing memory
    del dgdf
    
    # Calculating connected components
    _, conn_ids = connected_components(sparse_matrix)
    
    # Freeing memory
    del sparse_matrix 

    
    # Casting the two month GeoDataFrames together
    gdf = gpd.GeoDataFrame(pd.concat([gdf_prev_m, gdf_next_m], ignore_index = True),
                           crs = crs_m)
    
    # Freeing memory
    del gdf_prev_m, gdf_next_m
    
    # Erasing unneeded columns
    gdf = gdf.drop(columns = ['ID_prev', 'ID_next'])
    
    # Sanity check: length of gdf should be the same as that of conn_ids
    if len(conn_ids) != len(gdf):
        print('----------------------------------------------------------------')
        print('ERROR!!')
        print('The number of rows in prev_month-next_month dataframe is not \n the same as conn_ids!')
        print('')
        print('Exitting program.')
        print('----------------------------------------------------------------')
        sys.exit(1)
    
    
    # Otherwise, cast conn_ids as column identifying the connected components
    gdf['connected_ids'] = conn_ids
    
    # Freeing memory
    del conn_ids

    
    return gdf





def aggregate_connectedPolygons(gdf, num_process):
    """
    Function that casts into a MultiPolygon those polygons belonging to the 
    same group. Polygons belong to the same group if they share the same
    connected_ids value.
    
    This function works in parallel using dask-geopandas:
        https://dask-geopandas.readthedocs.io/en/latest/guide/dissolve.html
    accessed 25/07/2022

    Parameters
    ----------
    gdf : GeoDataFrame
        GeoDataFrame containing the polygons that belong to same conn_ids value 
        to be aggregated into MultiPolygons. The column conn_ids has repeated values.

    Returns
    -------
    gdf : GeoDataFrame
        GeoDataFrame containing the resulting MultiPolygons after casting together
        all those polygons that belong to the same group (have the same conn_ids
        value).

    """

    # First, create a dask-geopandas object from a geopandas one, with a certain
    # number of partitions
    ddf = dgpd.from_geopandas(gdf, npartitions=num_process)
    
    # Now, for the year 1993 we find a topology exception. I have identified the
    # polygon. I can either erase this polygon. But first I am going to add a buffer
    # to the geometries to see if I can circumvent the problem
    # ddf['geometry'] = ddf.geometry.buffer(0.5, join_style = 2).compute()
    
    # Prepare for dissolve the polygons sharing the same connected_ids value 
    # into MultiPolygons. Dissolve is like a groupby function. In this case,
    # we are grouping by 'connected_ids' and summarising to get the conn_ids, 
    # year, and month information. conn_ids and year are of course the same 
    # value for all polygons in group. But the months are two different oones. 
    # We keep the largest month because we will want to check if the resulting 
    # polygons are touching any of those from the following month.    
    shuffled = dissolve_shuffle(ddf, 'connected_ids', aggfunc={
        'conn_ids' : 'first',
        'year' : 'first',
         'month' : 'max',
     })
    
    # I believe the above function was just preparing the partitions and the operations
    # to carry out. Now we compute it.
    gdf = shuffled.compute()
    
    # COnverting again to a dask-geopandas dataframe
    # gdf = dgpd.from_geopandas(gdf, npartitions=num_process)
    
    # Then, undo the buffering
    # gdf['geometry'] = gdf.geometry.buffer(-0.5, join_style = 2)
    
    # The resulting polygons are in fact MultiPolygons of polygons within 150 m
    # but also of polygons that are touching one another, because the aobve 
    # procedure casts into a MultiPolygon, but does not dissolve the inner 
    # lines. We try to dissolve the inner lines
    # Getting rid of any left inner lines
    gdf['geometry'] = gdf['geometry'].buffer(0.001, join_style = 2)#.compute()
    # Finally, shrink back geometries again
    gdf['geometry'] = gdf['geometry'].buffer(-0.001, join_style = 2)#.compute()
    
    return gdf



def dissolve_shuffle(ddf, by=None, **kwargs):
    """
    Function that shuffle and partitions: I believe this means, shuffles the data
    so that elements sharing "by" (in this case , conn_ids value) are put together
    and then partitions between different "by" groups.
    
    This function was copied from 
    Taken from https://dask-geopandas.readthedocs.io/en/stable/guide/dissolve.html
    Under "Alternative solution"
    
    
    Parameters
    ----------
    ddf : Dask GeoDataFrame
        The GEoDataFrame that we want to dissolve.
    by : string, optional
        The column by which we want to group by and dissolve. The default is None.
    **kwargs : No idea
        DESCRIPTION.

    Returns
    -------
    I think it returns a map obkect linking the partitions and the function to
    be applied
        DESCRIPTION.

    """
    
    """Shuffle and map partition"""

    # Dissolve geometries within by group into a single geometry.
    meta = ddf._meta.dissolve(by=by, as_index=False, **kwargs)

    # I believe this bit organises the partitions to by conn_ids
    shuffled = ddf.shuffle(
        by, npartitions=ddf.npartitions, shuffle="tasks", ignore_index=True
    )

    # Finally, I believe this bit sets out that we want to dissolve the geometries
    # in parallel, with partitions organised around "by" (conn_ids)
    return shuffled.map_partitions(
        gpd.GeoDataFrame.dissolve, by=by, as_index=False, meta=meta, **kwargs
    )
        

"""
MAIN
"""
def main():
    
    # t.tic()
    
    # Handling warnings ------------------------------------------------------
    
    # Setting off slice DataFrame warning
    pd.options.mode.chained_assignment = None  # default='warn
    
    # Printing WARNING message that this script will erase the input files at the end
    # Asking if the user agrees
    print('-------------------------------------------------')
    print('WARNING: this script will erase the input files at the end:')
    print('../../data/temp/{}/MBfogo_c20_cerrado_polygons_{}_buff{}.{}')
    print('')
    user_decision = input('Are you OK with this? (y/n)')
    print('-------------------------------------------------')
    
    if user_decision == 'n':
        sys.exit(1)
    
    
    
    # User inputs ------------------------------------------------------------
    
    # Buffering distance. Needed to know what polygons to work on.
    # This buffer is not used in this script to merge polygons, only to indicate
    # which polygons to work on
    buffer_label = 150
    print('-------------------------------------------------')
    print('Using polygons from sameMonth obtained with a buffer of {} m'.format(buffer_label))
    print('-------------------------------------------------')
    
    # Asking the user whether to read server orders for a specific node
    server_order = input("Input '_node<>' to select a specific set of orders.\nFor no-selection, just press enter.\n")
    
    # CSV with the list of years to work with in this script
    csv_fp = '../server_orders/MBfogo_mergePolygons_consecutiveMonths_buffer{}.csv'.format(server_order) 
    
    # Reading list of files to work on (one after the other)
    list_years = fileList.readListFiles_fromCSV(csv_fp)
    
    # UPDATE 21/03/2023
    # I have observed that not all years have burned pixels in every month.
    # In particular, there are a few data files that do not contain data for 
    # the first few months. Hence, the start and end months are not set by the 
    # user, but by the data itself.
    # # The month we wish to start working on
    # start_month = 1
    # # The last month we want to work on
    # end_month = 12
    
    
    nprocess = 30
    
    # Merging polygons -------------------------------------------------------
    
    
    
    # Working on one file in server_order at a time (one year at a time)
    for year in list_years:
        
        print('')
        print('Working on year {}'.format(year))
        
        
        # ---------------- Input and output files ----------------------------
        # Input file with the polygons have been merged based on buffer and same-month
        # (after running them through fix_geometries)
        year_fp = '../../data/temp/{}/MBfogo_c20_cerrado_polygons_{}_buff{}.{}'
        # year_fp = '../../data/tests/beta/{}/MBfogo_c20_cerrado_polygons_{}_buff{}.{}'
        pol_fp = year_fp.format(year, year, buffer_label, 'shp')
        
        # Output directory
        out_dir = '../../data/temp/{}_CM/'.format(year)
        # out_dir = '../../data/tests/beta/{}_CM/'.format(year)
        # Creating output directory if it doesn't exist
        if not os.path.exists(out_dir.format(buffer_label)):
            os.makedirs(out_dir.format(buffer_label))
            
        
        # Output file
        shp_out_fp = out_dir + 'MBfogo_c20_cerrado_polygons_{}.shp'.format(year)
        
        # --------------------------------------------------------------------
        
        
        # Reading the full polygon's file with the data's information. 
        # This will be used just to get the index range containing a month's data
        # As well, there are years that do not contain any polygons for the first few months.
        # Hence, we also take the chance to get from the GeoDataFrame the information on the
        # start_month and end_month
        print('Reading polygons just to extract list of months')
        t.tic()
        month_indicator = gpd.read_file(pol_fp,
                                        include_fields = ['year','month'],
                                        ignore_geometry = True,
                                        engine = 'pyogrio')
        t.toc('Month indicator read in')
        # month_indicator = month_indicator[['year','month']]
        
        # Getting the months that are present in the data
        months_in_year = month_indicator['month'].unique()
        # print(months_in_year)
        # Getting the start and end months (the minimum and maximum)
        start_month = months_in_year.min()
        end_month = months_in_year.max()
        
        # Freeing memory
        # del months_in_year
        
        # Getting subset of polygons burned in start_month and the original data's CRS
        pols, crs_original = read_monthPolygons_convertCRS(pol_fp, month_indicator, start_month)

        
        
        # We work on one month's data at a time. Merging it with the previous month
        # Starting on second month, since start_month has no prior month
        # to merge it to
        for month in range(start_month + 1, end_month + 1):
            
            # It could be that there are polygons in month-1, but not in month
            # (e.g. data for the year 1995, there are no burnt pixels in February)
            if month not in months_in_year:
                # Then we save the month-1 data directly
                if not os.path.exists(shp_out_fp.format(buffer_label, year, buffer_label)):
                    # If the output file does not yet exist, create output file and save
                    pols.loc[pols['month'] == month - 1, ['conn_ids','year','month','geometry']].to_crs(crs_original).to_file(shp_out_fp.format(buffer_label, year, buffer_label))
                
                elif os.path.exists(shp_out_fp.format(buffer_label, year, buffer_label)):
                    # If file already exists, then append polygons to output file
                    pols.loc[pols['month'] == month - 1, ['conn_ids','year','month','geometry']].to_crs(crs_original).to_file(shp_out_fp.format(buffer_label, year, buffer_label), mode = 'a')

            elif month in months_in_year:
                
                # If month is in data, then we have to check if there are any
                # polygons in month-1. 
                # If there are no polygons in month-1, then we do not bother checking
                # month against month-1, convert pols as all month polygons and go on
                # to the next
                if month-1 not in months_in_year:
                    pols, _ = read_monthPolygons_convertCRS(pol_fp, month_indicator, month)
                    
                elif month-1 in months_in_year:
                    # If there is data for both months, then we have to check 
                    # the connected components and dissolve accordingly
            
                    print('')
                    print('Working on month {}'.format(month))
                    
                    # Reading this month's data and appending it to the previous month's
                    # polygons
                    
                    # Getting subset of polygons burned in month
                    t.tic()
                    pol_month, _ = read_monthPolygons_convertCRS(pol_fp, month_indicator, month)
                    t.toc('Polyons in month read in')
        
                    
                    # Appending polygons of this month (pol_month) to those of month-1 (pols)
                    # pols = gpd.GeoDataFrame(pd.concat([pols, pol_month], ignore_index = True), crs = 'epsg:5880')
                    
                    # Freeing memory
                    # del pol_month
                    
                    # pols carries over the polygons from the previous month
                    # pol_month contains the polygons of current month in loop
                    
                    # Identifying the connected components
                    print('Identifying connected polygons for months {} and {}...'.format(month, month - 1))
                    t.tic()
                    pols = id_connectedComponents(pols, pol_month, nprocess)
                    t.toc('Connected components identified in')
                    
                    # Freeing memory
                    del pol_month
        
                    
                    # Here, we should separate those polygons that are to be dissolved
                    # because they are touching one another, and those that are "isolated" 
                    # Regardin the "isolated" ones, those from month-1 can already be written
                    # to shapefile. We have to keep those "isolated" polygons of month
                    # because they could still be touching polygons in month+1.
                    print('Sorting non-connected components...')
        
                    # 1. Finding the duplicates (those polygons that are touching one 
                    #    another have been identified with the same connected_ids)
                    duplicates = pols['connected_ids'].duplicated(keep = False)
                    
                    
                    # 2. Hence, separating the "isolated" polygons from those polygons
                    # to be dissolved
                    pols_notDissolve = pols[~duplicates]
                    pols = pols[duplicates]
                    
                    t.tic()
                    # 3. Save and remove those polygons not to be dissolved from month-1
                    print('Saving non-connected polygons from previous month (month = {})'.format(month - 1))
                    if not os.path.exists(shp_out_fp.format(buffer_label, year, buffer_label)):
                        # If the output file does not yet exist, create output file and save
                        pols_notDissolve.loc[pols_notDissolve['month'] == month - 1, ['conn_ids','year','month','geometry']].to_crs(crs_original).to_file(shp_out_fp.format(buffer_label, year, buffer_label))
                    
                    elif os.path.exists(shp_out_fp.format(buffer_label, year, buffer_label)):
                        # If file already exists, then append polygons to output file
                        pols_notDissolve.loc[pols_notDissolve['month'] == month - 1, ['conn_ids','year','month','geometry']].to_crs(crs_original).to_file(shp_out_fp.format(buffer_label, year, buffer_label), mode = 'a')
                    
                    # Now we can removing these polygons (so, we keep those only those that are in month)
                    pols_notDissolve = pols_notDissolve[pols_notDissolve['month'] == month]
                    
                    t.toc('Non-connected polygons dealt with in')
        
                    
                    # 4. Dissolve connected polygons from month and month-1 to a MultiPolygon
        
                    print('Dissolving polygons from month {} into MultiPolygons...'.format(month))
                    t.tic()
                    pols = aggregate_connectedPolygons(pols, nprocess)
                    t.toc('Connected polygons aggreated in')
                    
                    # 5. Finally, concatenating back these two set of polygons to see 
                    # if they are connected to any next month's polygons
                    pols = gpd.GeoDataFrame(pd.concat([pols_notDissolve, pols], ignore_index = True), crs = 'epsg:5880')
                    pols = pols[['conn_ids','year','month','geometry']]
                    
                    # Freeing memory
                    del pols_notDissolve
            
            
            print('')
            
        
        print('Saving last polygons to file')
        # Saving the polygons left from end_month, appending them to file
        pols.to_crs(crs_original).to_file(shp_out_fp.format(buffer_label, year, buffer_label), mode = 'a')
           
        
        
        print('\nErasing input directory...')
        in_dir = os.path.dirname(year_fp).format(year)
        # print(in_dir)
        try:
            shutil.rmtree(in_dir)
        except OSError as e:
            # If it fails, inform the user.
            print("Error: %s - %s." % (e.filename, e.strerror))

        print('---------------------------------------------------------------')
        print('\n\n')
        
    t.toc('Programme finished in')    
    
        
    # Setting on slice DataFrame warning
    pd.options.mode.chained_assignment = 'warn'

if __name__ == '__main__':
    main()
    
