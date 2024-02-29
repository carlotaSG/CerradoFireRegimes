# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 10:44:47 2021

@author: scat8298
"""
import sys
import os
from multiprocessing import Pool


def call_fix(input_fp, output_fp):
    with open('./fix_geometries.py') as fix:
        
        script = fix.read()
        
        print(os.getcwd())
        
        sys.argv = ['./fix_geometries.py',  input_fp, output_fp]
        
        exec(script)

def main():
    
    poolInput = [('../../data/MBfogo_c10_cID021290_2019.shp', '../../data/fixing_server_cID021290.shp'),
                 ('../../data/MBfogo_c10_cID021289_2019.shp', '../../data/fixing_server_cID021289.shp')]
    
    
    with Pool(2) as p:

        print('Starting parallel computation within calling_fix.')
        
        # Sending to polygonise in parallel
        # Return list of tuples: GeoDataFrames per cell and cell identifiers
        p.starmap(call_fix, poolInput)
    
    
if __name__ == '__main__':
    main()