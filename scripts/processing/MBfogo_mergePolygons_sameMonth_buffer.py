# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 10:21:32 2022

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that joins spatially, into a single MultiPolygon, those MBfogo polygons 
that belong to the same month and that are within a certain distance (buffer).
The files (years) to work on are specified by the user in a server order file.

This script carries out the task working over each month sequentially and 
using dask-geopandas to parallelise over the polygons to be dissolved within 
each month.

WARNING: this script should be run in the server because, to save disk space, 
         the input files are erased at the end.

Input files:
    ../../data/processed/MBfogo_c{}_cerrado_polygons/MBfogo_c{}_cerrado_polygons_<year>.shp
    ../server_orders/MBfogo_mergePolygons_sameMonth_buffer<buffer_meters>.csv
    
Output files:
    ../../data/temp/<year>/MBfogo_c{}_cerrado_polygons_<year>_buff<buffer_meters>.shp
    
    
This script contians the following functions:
    * id_connectedComponents - identifies with the same ide value the groups of 
                               polygons that are within a buffer distance from one another
    * aggregate_connectedPolygons - casts into a MultiPolygon those pixels sharing 
                                    the same conn_ids value. Uses dask-geopandas
                                    and it calls dissolvve_shuffle
    * dissolve_shuffle - prepares the data into partitions for polygons to be 
                         cast into MultiPolygons using dask-geopandas (this is 
                         done in aggregate_connectedPolygons)
"""

"""
Imports
"""
import geopandas as gpd
import pandas as pd
import dask_geopandas as dgpd
from libpysal.weights import fuzzy_contiguity
import os, glob


from pytictoc import TicToc
t = TicToc()

import sys
sys.path.append('../utils/')

import fileList


"""
Functions
"""

def id_connectedComponents(gdf, buff, last_W):
    """
    Function that receives a GeoDataFrame and identifies those polygons that 
    are within a certain buffer distance from one another.

    Parameters
    ----------
    gdf : GeoDataFrame
        The GeoDataFrame containing the polygons of which to identify the connected
        components.
    buff : int
        Distance to be used to consider polygons as being connected.
    last_W : int
        At some point I decided that the output of the full script would be a 
        GeoDataFrame where each polygon would be identified by a single ID.
        Because the script works sequentially over months, the id of the polygons
        needs to be carried over from the previous month. Variable last_W is the 
        last/largest id value from the previous month.

    Returns
    -------
    gdf : GeoDataFrame
        The input GeoDataFrame with an additional column 'conn_ids' that identifies
        all groups of polygons found within a buff distance with the same number.
    last_W : int
        The new largest id assigned to a group of polygons. This id will be 
        inherited by this same function to update the poylgons ids.

    """

    # Identifying the polygons that are within a distance from one another
    W = fuzzy_contiguity(gdf, buffering = True, buffer = buff, silence_warnings = True)
        
    # Correcting to accumulate index through months
    conn_ids = W.component_labels + last_W + 1
    
    # Casting this id vector as a column
    gdf['conn_ids'] = conn_ids
        
        
    # Updating the last weight (polygon id) position
    last_W = max(conn_ids)
        
    # Freeing memory
    del W, conn_ids
        
    
    return gdf, last_W





def aggregate_connectedPolygons(gdf):
    """
    Function that casts into a MultiPolygon those polygons belonging to the 
    same group. Polygons belong to the same group if they share the same
    conn_ids value.
    
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
    ddf = dgpd.from_geopandas(gdf, npartitions=38)
    
    # Dissolve the polygons sharing the same conn_ids value into MultiPolygons
    shuffled = dissolve_shuffle(ddf, 'conn_ids')
    
    # I believe the above function was just preparing the partitions and the operations
    # to carry out. Now we compute it.
    gdf = shuffled.compute()
      
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
    t.tic()
    
    # The MapBiomas collection we are working with
    collection = 2
    
    # Setting off slice DataFrame warning
    pd.options.mode.chained_assignment = None  # default='warn
    
    # Printing WARNING message that this script shall be ran only in the server
    # because this script erases original files to free memory
    print('-------------------------------------------------')
    print('THIS SCRIPT SHALL ONLY BE RAN IN THE LINUX CLUSTER!')
    print('')
    print('Because it erases the original files:')
    print('../../data/processed/MBfogo_c{}0_cerrado_polygons/MBfogo_c{}0_cerrado_polygons_{}.{}'.format(collection, collection, '{}', '{}'))
    print('')
    user_decision = input('Are you OK with this? (y/n)\n')
    print('-------------------------------------------------')
    
    if user_decision == 'n':
        sys.exit(1)
    
    # User inputs ------------------------------------------------------------
    
    # Select the buffering distance to consider two polygons as being part of the
    # same MultiPolygon (meters)
    
    buffer_label = 150
    buffer = int(buffer_label/2)
    
    print('')
    print('-------------------------------------------------')
    print('Using a buffer of {} m'.format(buffer_label))
    print('Imposing {} as a one-sided buffer.'.format(buffer))
    print('-------------------------------------------------')
    print('')
    
    
    # Asking the user whether to read server orders for a specific node
    server_order = input("Input '_node<>' to select a specific set of orders.\nFor no-selection, just press enter.\n")
    
    # CSV with the list of raster files to work with in this script
    csv_fp = '../server_orders/MBfogo_mergePolygons_sameMonth_buffer{}.csv'.format(server_order) 
    
    # Reading list of files to work on (one after the other)
    list_years = fileList.readListFiles_fromCSV(csv_fp)
    
    
    # Working on one file in server_order at a time (one year at a time)
    for year in list_years:
        
        print('')
        print('Working on year {}'.format(year))
        
        
        # Input directory with original files
        year_fp = '../../data/processed/MBfogo_c{}0_cerrado_polygons/MBfogo_c{}0_cerrado_polygons_{}.{}'.format(collection, collection, '{}', '{}')
        pol_fp = year_fp.format(year, 'shp')
    
        
        
        # Output directory
        out_dir = '../../data/temp/{}/'.format(year)
        # Output file
        shp_out_fp = out_dir + 'MBfogo_c{}0_cerrado_polygons_{}_buff{}.shp'.format(collection, year, buffer_label)
        
        
        # Reading the shapefile to get the index-month relation to later be 
        # able to retrieve only those polygons mapped in a certain month
        print('Reading polygons just to extract list of months')
        
        # UPDATE 11/10/2023 using csv file to accelerate SHP reading
        data_csv = gpd.read_file(pol_fp)#, rows = 1000)
        data_csv = data_csv[['year','month']]
        
        
        # Initialising id counter
        # This counter will allow each polygon in the output to have a unique id
        # as it will hep keep track of the last id (weight) value from one month
        # to the next
        last_W = -1
        
        
        # If output folder does not exist, create it
        if not os.path.exists(out_dir.format(buffer_label)):
            os.makedirs(out_dir.format(buffer_label))
        
        
        # Working on a month by month basis
        for month in range(1, 12 + 1):
            
            print('Working on month {}'.format(month))
            
            # Get the row indices that contain this month's data
            month_index = data_csv[data_csv['month'] == month].index
            
            # It could be that there are no pixels burnt in this month.
            # In that case we just jump to the next year
            if len(month_index) == 0:
                print('WARNING: no pixels burnt in month {} of year {}'.format(month, year))
            elif len(month_index) > 0:
            
                # Reading and subsetting relevant data
                pol = gpd.read_file(pol_fp, rows = slice(month_index[0], month_index[-1]))
                
                # Freeing memory 
                del month_index
    
    
                # Converting to metre CRS
                print("Converting polygon's CRS")
                crs_original = pol.crs
                pol = pol.to_crs('epsg:5880')
            
            
                # Identifying the groups of polygons that are within a buffering
                # distance from one another
                print('Identifying connected polygons for month {}...'.format(month))
                connected_pol, last_W = id_connectedComponents(pol, buffer, last_W)
                
                # Freeing memory
                del pol
            
                # Here, we should separate those polygons found in a group of polygons 
                # that are within a buffering distance from one another - and hence 
                # that have to be dissolved -, from those that are isolated - and so 
                # they can be written to output file directly. 
                print('Sorting non-connected components...')
    
                # 1. Finding the duplicates - as a boolean vector indicating those
                #    polygons belonging to a group of connected polygons (True)
                #    from those that are isolated (False)
                duplicates = connected_pol['conn_ids'].duplicated(keep = False)
                
                
                # 2. Hence, separating those polygons not to be dissolved, and those that are
                pols_notDissolve = connected_pol[~duplicates]
                connected_pol = connected_pol[duplicates]
                
                
                print('Saving non-connected polygons for month {}...'.format(month))
                # Further, save and remove those polygons not to be dissolved from month-1
                if not os.path.exists(shp_out_fp.format(buffer_label, year, buffer_label)):#month == 1:
                    # If we are working on the first round of month pairs, create output file and save
                    pols_notDissolve[['conn_ids','year','month','geometry']].to_crs(crs_original).to_file(shp_out_fp.format(buffer_label, year, buffer_label))
                elif os.path.exists(shp_out_fp.format(buffer_label, year, buffer_label)): #else:
                    # Write these polygons appending them to output file
                    pols_notDissolve[['conn_ids','year','month','geometry']].to_crs(crs_original).to_file(shp_out_fp.format(buffer_label, year, buffer_label), mode = 'a')
                
                # Freeing memory
                del pols_notDissolve
    
            
                # If there are any groups of polygons left to connect, aggregate to 
                # a MultiPolygon
                if len(connected_pol) > 0:
                    print('Dissolving polygons into MultiPolygons...')
                    connected_pol = aggregate_connectedPolygons(connected_pol)
    
                
                    # Now, save the shapefile with the connected polygons
    
                    print('Saving files for month {}'.format(month))
                    connected_pol = connected_pol.to_crs(crs_original)
                    
                    if not os.path.exists(shp_out_fp.format(buffer_label, year, buffer_label)):
                        connected_pol[['conn_ids','year','month','geometry']].to_file(shp_out_fp.format(buffer_label, year, buffer_label))
                    elif os.path.exists(shp_out_fp.format(buffer_label, year, buffer_label)):
                        connected_pol[['conn_ids','year','month','geometry']].to_file(shp_out_fp.format(buffer_label, year, buffer_label), mode = 'a')
    
                
                del connected_pol
            
            print('')
            
        
        print('\nRemoving input files for year {}'.format(year))
        # Use library glob to return all files that match the pathname with the wildcard
        for file in glob.glob(year_fp.format(year, '*')):
            print(file)
            try:
                os.remove(file)
            except OSError as e:
                # If it fails, inform the user.
                print("Error: %s - %s." % (e.filename, e.strerror))
            
        
        print('Finished merging polygons of the same month with spatial criteria.')
        
        print('\n \n \n ----------------------------------------------')
        
        
    t.toc('Program ran in')
        
    # Setting on slice DataFrame warning
    pd.options.mode.chained_assignment = 'warn'

if __name__ == '__main__':
    main()

