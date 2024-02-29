# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 18:01:54 2023

@author: scat8298
"""
import pandas as pd

for year in range(1985, 2020 + 1):
    data = pd.read_csv('../../data/processed/climate_data/terraclimate_grid50km_{}.csv'.format(year))
    
    data['year'] = year
    
    data = data[['id','year','month','mean_prec','mean_aet','mean_pet']]
    
    data.to_csv('../../data/processed/climate_data/terraclimate_grid50km_{}.csv'.format(year), index = False)