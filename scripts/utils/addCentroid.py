# -*- coding: utf-8 -*-
"""
Created on Fri May 27 12:08:10 2022

@author: scat8298

Add ceontroid coordinates
"""

import geopandas as gpd

gpd_fp = '../data/grids/cerrado_grid_02deg_cropped.shp'

grid = gpd.read_file(gpd_fp)


grid["lon"] = grid.centroid.x
grid["lat"] = grid.centroid.y

grid.to_file(gpd_fp)