# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 14:00:13 2022

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculates the fraction of area occupied by polygons overlapping a
geometry. It finds those polygons overlapping the geometry, clips them, and 
calculates the area fraction of the clipped polygons.

It can be run for a list of polygons and a unique geometry, or for a list of 
polygons and geometries. In the latter case, it calculates the fraction of areas
for each geometry in list.

"""


"""
Imports
"""
import pandas as pd
import geopandas as gpd
from multiprocessing import Pool

"""
Functions
"""
# def to_m_CRS():
    
def clip_calculateArea(gdf_polygons, geom, pol_id, geom_id):
    
    # where pol_id is the column of gdf_polygons with the ID for the polygons
    # where geom_id is the column of geom with the ID of the geom
    # returns a dataframe: geom_id | pol_id | frac_areas
    
    # # Checking that both GeoDataFrames are in teh same CRS, and that its units
    # # are meters
    # if (gdf_polygons.crs == gdf_geometry.crs) & (gdf_polygons.crs['units'] == )
    # geom = geom.to_crs(gdf.crs)
    
    # Setting off warning momentarily
    pd.options.mode.chained_assignment = None  # default='warn'
    
    # Calculating area of geometry
    geom_area = geom.area
    
    # Obtaining the list of selected polygons and their clipped geometries
    list_clipped_gdf = gpd.clip(gdf_polygons['geometry'], geom)
    
    # GeoDataFrame of selected polygons
    clipped_gdf = gdf_polygons[gdf_polygons.index.isin(list_clipped_gdf.index)]
    
    # Inserting the clipped geometries
    clipped_gdf['geometry'] = list_clipped_gdf

    # Calculating areas of clipped polygons
    clipped_gdf['area_km2'] = clipped_gdf.geometry.area
    
    # Calculating fraction occupied by each polygon
    clipped_gdf['frac_area'] = clipped_gdf['area_km2'] / geom_area
    
    
    # Formatting output
    clipped_gdf = pd.DataFrame({'geom' : geom_id, 
                                pol_id : clipped_gdf[pol_id], 
                                'frac_area' : clipped_gdf['frac_area']})
    
    # Setting on warning again
    pd.options.mode.chained_assignment = 'warn'  # default='warn'
    
    return clipped_gdf
    

def clip_calculateArea_withPool(gdf_polygons, gdf_geom, pol_id, geom_id, num_process):
    # Runs the function clip_calculateArea parallelizing over the geometries in gdf_geom
    # pol_id indicates the column of gdf_polygons containing the polygons' ID
    # geom_id indicates the column of gdf_geom containing the geometries' ID
    # print(gdf_geom)
    # print([i for i in range(len(gdf_geom))])
    
    
    
    # Prepare input to polygonise in parallel using map
    poolInput = [(gdf_polygons, gdf_geom.loc[i, 'geometry'], 
                  pol_id, gdf_geom.loc[i, geom_id]) for i in range(len(gdf_geom.geometry))]

    
    # polygonise raster by parts (cells in grid) and working in parallel
    with Pool(num_process) as p:

        print('Starting parallel computation within clip_calculateArea_withPool.')
        
        # Sending to polygonise in parallel
        # Return list of tuples: GeoDataFrames per cell and cell identifiers
        output = p.starmap(clip_calculateArea, poolInput)
        
    # Freeing some memory (erasing input to parallelisation)
    del poolInput
    
    # Formatting output
    output = pd.concat(output, ignore_index = True)
    # print(output)
    
    return output
    
    
    
    
    
    
    
    
    
    


