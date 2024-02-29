# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 17:05:36 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that merges raster files cooresponding to different parts of the same area
into a single raster file. It does so for two different types of MapBiomas [FOGO]
data: burned area, and frequencies.

This script must be run from the command line along with one argument indicating
which program option is to be run:
    - Option 1: merge MapBiomas [FOGO] burned area files. Data must be stored in 
    '../../data/MBfogo_c10/' in different folders depending on the year. Output
    written to '../data/raw/MBfogo_ba/'
    
    - Option 2: merge MapBiomas [FOGO] frequency data. Data must be stored in 
    '../../data/MBfogo_frequencies_1985-2020/'. Output goes to '../data/raw/MBfogo_freq/'
"""

"""
Imports
"""
import os
import rasterio

import sys
sys.path.append('../utils/')
import fileList, merge_rasters


"""
Functions
"""
def merge_MBfogo_ba():
    """
    If option 1 is chosen, this function runs. It reads the folder structure
    from the input_fp, which should be a list of directories with a year as directory name.
    Then, it loops over the year directories reading the list of .tif files inside,
    merges them into a single raster file, and saves them to output_fp with a
    respective filename.

    Returns
    -------
    None.

    """
    
    print("Merging MapBiomas [FOGO] burned area data.")
    
    # input_fp = '../../../data/MBfogo_c10/'
    # input_fp = 'D:/carlota/projects/data/MBfogo_c20'
    input_fp = 'C:/Users/scat8298/OneDrive - Nexus365/Documents/carlota/projects/_chapter02/data/raw/MBfogo_ba'
    

    list_dir_fp = os.walk(input_fp)
    
    # List of directories each containing the raster subfiles for one year in the
    # MapBiomas dataset
    list_dir_fp = [element.path for element in os.scandir(input_fp) if element.is_dir()]
    
    # For each directory, get list of child filepaths and merge into one raster, then save
    for dir_fp in list_dir_fp:
        print('Working on year {}'.format(os.path.basename(dir_fp)))
        
        list_raster_fp = fileList.fileList(dir_fp)
        
        out_raster, out_metadata = merge_rasters.merge_raster(list_raster_fp)
        
        # Save as raster
        # First, create output filepath
        year = os.path.basename(dir_fp)
        # output_fp = '../../data/raw/MBfogo_ba/MBfogo_c10_cerrado_{}.tif'.format(year)
        # output_fp = 'D:/carlota/projects/data/MBfogo_c20/MBfogo_c20_cerrado_{}.tif'.format(year)
        output_fp = 'C:/Users/scat8298/OneDrive - Nexus365/Documents/carlota/projects/_chapter02/data/raw/MBfogo_ba/MBfogo_c20_cerrado_{}.tif'.format(year)

        print('Writing tif file.')
        # Save tif file
        with rasterio.open(output_fp, 'w', **out_metadata) as output:
            output.write(out_raster)
        
        print('\n')
        #break
        
def merge_MBfogo_freq():
    """
    If option 2 is chosen, this function runs. It reads the list of tif files 
    within the input_fp folder, merges them into a single raster file, and 
    saves them to output_fp with a respective filename.

    Returns
    -------
    None.

    """
    
    # TODO: make this more flexible so that different frequency files could be
    # merged consecutively
    print('Merging fire frequencies for period 1985 - 2020')
    
    # Input folder
    input_fp = '../../../data/MBfogo_frequencies_1985-2020/'
    
    # Reading list of files in folder
    list_fp = fileList.fileList(input_fp)
    
    # Merging raster files into one
    out_raster, out_metadata = merge_rasters.merge_raster(list_fp)
    
    # Output filepath    
    output_fp = '../../data/raw/MBfogo_freq/MBfogo_freq_c10_cerrado_1985-2020.tif'
    
    # Saving file
    with rasterio.open(output_fp, 'w', **out_metadata) as output:
        output.write(out_raster)
        
    print('Finished merging MapBiomas [FOGO] frequency raster files.')
    
    
    

"""
MAIN
"""
def main():
    
    # # User command line input chooses the mode of script to run:
    # # Option 1: Input == 1, merge MapBiomas [FOGO] burned area files 
    # # Option 2: Input == 2, merge MapBiomas [FOGO] frequency files
    
    # # User must give one input argument
    # try:
    #     prog_mod = int(sys.argv[1])
    # except:
    #     print('ERROR')
    #     print("Program option required as input. None introduced.")
    #     sys.exit(1)
        
    
    # # Depending on user input, a program option is run
    # if prog_mod == 1:
    #     merge_MBfogo_ba()
        
    # elif prog_mod == 2:
    #     merge_MBfogo_freq()
    
    merge_MBfogo_ba()
        
    # If user input is netiher 1 nor 2, we notify of an error and finish program normally.
    # else:
    #     print('ERROR.')
    #     print('Wrong program option input. Options avaialable are 1 for burned area, and 2 for frequencies.')
    
    
if __name__ == '__main__':
    main()