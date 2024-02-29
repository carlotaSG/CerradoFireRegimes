# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 20:00:01 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that reads all files that contain data per cell for different periods. 
Returns all data in a single file, along with the geometry information,
having filtered out cells with less than a certain percentage of their area inside the Cerrado (optionally).
"""

"""
Imports
"""
import pandas as pd
from add_gridInfo import add_gridInfo

"""
Functions
"""


def readFormat_allPeriods(data_fp, nperiods, grid_id, thresh = None, incl_cells = None, excl_cells = None):
    
    
    df = []

    for i in range(1, nperiods + 1):
        df_i = pd.read_csv(data_fp.format(i))
        df_i = add_gridInfo(df_i, grid_id, filter_thresh = thresh, incl_cells = incl_cells, excl_cells = excl_cells)
        df.append(df_i)
        del df_i

    df = pd.concat(df, ignore_index = True)
    
    return(df)
