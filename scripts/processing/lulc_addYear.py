# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 18:01:54 2023

@author: scat8298
"""
import pandas as pd

grid_id = '30km'

for year in range(1985, 2020 + 1):
    data = pd.read_csv('../../data_new_grids/processed/lulc_data/MBlanduse_percArea_subsetLULC_grid{}_{}.csv'.format(grid_id, year))
    
    data['year'] = year
    
    data = data[['id','year','anthropic','natural','flammable','nonflammable', '4', '3', '12', '11', '15', '18','9','21']]
    
    data.to_csv('../../data_new_grids/processed/lulc_data/MBlanduse_percArea_subsetLULC_grid{}_{}.csv'.format(grid_id, year), index = False)