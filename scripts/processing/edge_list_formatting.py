# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 16:34:52 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that gets a distance matrix between polygons (or other geometric elements)
and formats it to be used as an "edge list" input file for spatially constrained
hierarchical clustering. It include only those edges whose distance is smaller than
a certain threshold. It discard the edges of distance 0 (a polygons with itself).


Input format:
    id1 | id2 | dist_m

Output format:
    from | to | distance
    

"""

"""
Imports
"""
import pandas as pd


"""
MAIN
"""
if __name__ == '__main__':
    
    # ------------------- User inputs ----------------------------------------
    
    
    # Grid identifier
    grid_id = '50km'
    
    # Input distance matrix filepath
    dist_mtx_fp = '../../data/processed/clustering_algorithm/distance_matrix_grid_{}_epsg5880.csv'.format(grid_id)
    # dist_mtx_fp = '../../data/processed/clustering_algorithm/distance_matrix_80cells_grid_{}_epsg5880.csv'.format(grid_id)
    
    # Maximum distance of separation between polygons below which we select the
    # polygons pairs
    if grid_id == '50km':
        # TIP: suggesting adding a few meters to account for inaccuracies in the distances
        dist_thresh = 70710 + 100
    elif grid_id == '30km':
        dist_thresh = 42426 + 100
    # This distance implies a certain connectivity (the number of cells in grid considered
    # as being consecutive)
    connectivity = 8
    
    # Output condensed distance matrix
    out_dst_mtx_fp = '../../data/processed/clustering_algorithm/edge_list_connectivity{}_grid_{}_epsg5880.csv'.format(connectivity, grid_id)
    # out_dst_mtx_fp = '../../data/processed/clustering_algorithm/edge_list_connectivity{}_80cells_grid_{}_epsg5880.csv'.format(connectivity, grid_id)
    
    
    # ------------- Formatting distance matrix -------------------------------
    
    # Read original distance matrix
    dist_mtx = pd.read_csv(dist_mtx_fp)
    
    # Subsetting only distances smaller than thresh and non-zero
    dist_mtx = dist_mtx[(dist_mtx['dist_m'] > 0) & (dist_mtx['dist_m'] < dist_thresh)]
    
    # Formatting the name of the columns
    dist_mtx.columns = ['from', 'to', 'distance']
    
    
    # ----------------- Saving output ----------------------------------------
    
    dist_mtx.to_csv(out_dst_mtx_fp, index = False)
    
    
    