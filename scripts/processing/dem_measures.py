# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 09:32:50 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculates the slope, roughness and Terrain Ruggedness Index from 
an SRTM DEM file.


"""

"""
Imports
"""
import os
from osgeo import gdal
from pytictoc import TicToc
t = TicToc()


"""
MAIN
"""

def main():
    
    # ------------------- User inputs ----------------------------------
    
    measures = ['roughness', 'TRI', 'slope']
    
    input_file = '../../data/raw/SRTM_DEM/SRTM_DEM_merged.tif'
    
    output_file = '../../data/processed/SRTM_DEM/SRTM_DEM_merged_{}.tif'
    
    cmd = "gdaldem {} {} {} -compute_edges"
    
    # ------------ Processing ------------------------------------------
    
    # Working on a topographic measure at a time
    for measure in measures:
        
        print('Working on topographic measure: {}'.format(measure))
        
        measure_out_file = output_file.format(measure)
        
        measure_cmd = cmd.format(measure, input_file, measure_out_file)
        
        # Calculate the output file
        t.tic()
        os.system(measure_cmd)
        t.toc('{} produced in'.format(measure))
        
        print('')
    

    
    
    
    
if __name__ == '__main__':
    main()
    

