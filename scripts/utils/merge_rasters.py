# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 17:49:55 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that receives a list of raster filepaths that merges one by one and returns
the raster and the transform
"""

"""
Imports
"""
import rasterio.merge


"""
Functions
"""
def merge_raster(list_raster_fp):
    
    print('Reading raster files...')
    
    # List where to store the different raster datasets opened in 'r' mode
    list_raster = []
    for raster_fp in list_raster_fp:
        #print('Opening raster file: {}'.format(raster_fp))
        
        # Open in 'r' mode
        raster_dataset = rasterio.open(raster_fp)
        
        
        
        # Append new raster to list
        list_raster.append(raster_dataset)
        
    print('Merging files...')
    # Merge raster datasets
    merged_raster, merged_transform = rasterio.merge.merge(list_raster)
    
    # Update metadata:
    # First, get original metadata
    metadata = raster_dataset.meta
    # Update new fields
    metadata.update({
        'width': merged_raster.shape[2],
        'height': merged_raster.shape[1],
        'transform': merged_transform
        })
    
    
    # Close raster datasets to avoid further problems
    for raster_dataset in list_raster:
        raster_dataset.close()
    
    
    return merged_raster, metadata
