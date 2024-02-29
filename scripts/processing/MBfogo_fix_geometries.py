# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 11:48:32 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that reads in a list of shapefiles, fixes the geometries using PyQGIS and
ovwerwrites the files.

Input directories and output directories are indicated in a server_orders file: first column
for input files, second column for output files. This server_orders file must be 
updated depending on the job that needs doing. Updating is done manually, as 
server_order file is hard-coded in the script.

Script can run sequentially (one file after another) or in parallel. But running
in parallel is just correcting invalid geometries parallelising over the SHP files
found in the same input folder. That is, it does not parallelise over the polygons
within the SHP file. This is because the PyQGIS script called by this script 
(./scripts/utils/fix_gemoetries.py) operates over the full file.


WARNING: input files are deleted at the end!

Function in this script:
    * call_fix - prepares the input and executres the script ../utils/fix_geometries.py

"""


"""
Imports
"""
import os
import shutil
from multiprocessing import Pool

import sys
sys.path.append('../utils')

import fileList
# import fix_geometries

from pytictoc import TicToc
t = TicToc()

"""
Functions
"""
def call_fix(input_fp, output_fp):
    """
    Function that prepares the input and executes the script ../utils/fix_geometries.py,
    which uses PyQGIS to check polygon's validity and corrects them.

    Parameters
    ----------
    input_fp : string
        Filepath to shapefile of geometries that we want to fix.
    output_fp : string
        Filepath to shapefile to save fixed geometries.

    Returns
    -------
    None.

    """
    with open('../utils/fix_geometries.py') as fix:
        
        # Read script
        script = fix.read()
        
        # Setting arguments to pass to fix_geometries execution
        sys.argv = ['../utils/fix_geometries.py',  input_fp, output_fp]
        
        # Call execution
        exec(script)


"""
MAIN
"""

def main():
    t.tic()
    
    # Printing WARNING message that this script will erase the 
    # input files and directory at the end because these are no longer needed.
    print('-------------------------------------------------')
    print('WARNING!!')
    print('')
    print('This script erases the input files and directories at the end.')
    print('')
    user_decision = input('Are you OK with this? (y/n)\n')
    print('-------------------------------------------------')
    
    if user_decision == 'n':
        sys.exit(1)
    
    #---------------------- User inputs -------------------------------
    # Reading list of directories over which to operate - one at a time
    # Asking the user whether to read server orders for a specific node
    server_order = input("Input '_node<>' to select a specific set of orders.\nFor no-selection, just press enter.\n")
    
    # CSV with the list of directories (one per year) to work with in this script
    # csv_fp = '../server_orders/MBfogo_fix_geometries_after_consecutiveMonthsBuffer{}.csv'.format(server_order)
    csv_fp = '../server_orders/MBfogo_fix_geometries{}.csv'.format(server_order)
    
    
    
    # Decide on whether to work sequentially (s) or in parallel (p)
    mode = 'p' #p s
    # Number of workers
    num_process = 10
    
    # UPDATE 16/03/2023: Output directory is now indicated in the server_orders
    # file
    # Output directory
    # HINT: The output filepath where to write the fixed geometries MUST BE 
    # different from the input. Hence, we use a temporary folder that we create,
    # and the same filepath (see below function fileList)
    # out_dir = '../../data/temp/fixed/{}' # the {} is for the year, set in the loop
    
    #---------------------- The programme -----------------------------
    
    # Reading list of raster filepaths
    list_dir = fileList.readListFiles_fromCSV(csv_fp)
    
    print(list_dir)
    
    # for (in_dir, out_dir) in list_dir:
    for in_dir in list_dir:
        
        out_dir = '../../data/processed/MBfogo_c20_cerrado_polygons_b150_CM_fixed'
    
        # Input directory
        # in_dir = '../data/polygons_tocheck/2019'
        
        # UPDATE 16/03/2021: Now we print directly the name of the input directory
        # The year that we are working in is given by the name of the folder 
        # we are working on
        # year = _ath.basename(in_dir)
        
        print('-----------------------------------------------------------')
        print('Working on {}'.format(in_dir))

        # UPDATE 16/03/2023
        # # Output directory
        # out_dir = out_dir.format(year)

        # If the desired output folder does not exist, create it
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        
        # Read list of files for which to fix the geometries, generate output filepaths 
        # to be the same as input.
        list_fp = fileList.fileList(in_dir, outDir = out_dir, in_extensions = 'shp')
        list_fp = list(zip(list_fp[0], list_fp[1]))

        # Sequentially 
        if mode == 's':
            # One file at a time
            for in_fp, out_fp in list_fp:
                
                print('Working on file {}'.format(os.path.basename(in_fp)))
                t.tic()
                call_fix(in_fp, out_fp)
                t.toc('Geometries fixed in')
                
                
        # In parallel
        if mode == 'p':
            with Pool(num_process) as p:
            
                print('Starting parallel computation within calling_fix.')
                
                t.tic()
                p.starmap(call_fix, list_fp)
                t.toc('Geometries fixed in')
    
        
        # Finally, copy all those .cpg files that are not rewritten by QGIS when
        # fixing geometries to the desired folder
        # First, select /cpg files
        cpg_files = fileList.fileList(in_dir, outDir = out_dir, in_extensions = 'cpg')
        cpg_files = list(zip(cpg_files[0], cpg_files[1]))
        # Copy them one by one
        for in_file, out_file in cpg_files:
            shutil.copy(in_file, out_file)
            
        # UPDATE 16/03/2023
        # Finally, we erase the input directory and files as we are no longer
        # interested in having files with potenitally invalid geometries
        print('\nErasing input directory...')
        try:
            shutil.rmtree(in_dir)
        except OSError as e:
            # If it fails, inform the user.
            print("Error: %s - %s." % (e.filename, e.strerror))
        
        print('\n\n')

    
if __name__ == '__main__':
    main()

