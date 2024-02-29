# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 11:07:32 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that performs data imputation to those cells in a period with less than X
fire polygons for a certain list of variables. It read the cell data for a certain 
period to identify the cells with <= X polygons, identifies the neighbours using the
edge list, performs imputation by assigning to the cell the average value of those
8 neighbours with more than X polygons. 

Input files:
    ../../data/processed/clustering_algorithm/input_vars_grid{}_period{}.csv
    '../../data/processed/clustering_algorithm/edge_list_connectivity8_grid_{}_epsg5880.csv'
    
Output files:
    ../../data/processed/clustering_algorithm/input_vars_imputation{}_grid{}_period{}.csv
    
"""

"""
Imports
"""
import pandas as pd

"""
Functions
"""
def get_neighbours(cell, data, connections, thresh):
    
    # First, find neighbours
    candidate_neighbours = connections.loc[(connections['from'] == cell) | (connections['to'] == cell), ['from', 'to']]
    candidate_neighbours = candidate_neighbours['from'].tolist() + candidate_neighbours['to'].tolist()
    candidate_neighbours = [i for i in candidate_neighbours if i != cell]
    
    # Subset data for the candidate neighbours if tey have enough number of fires
    data_subs = data[(data['id'].isin(candidate_neighbours)) & (data['number_fires'] > thresh)]
    
    # Final list of neighbours
    candidate_neighbours = data_subs['id'].tolist()
    
    return(candidate_neighbours)


"""
MAIN
"""

def main():
    
    # -------------------------- User inputs ----------------------------------
    
    # Grid id
    grid = '50km'
    
    # Periods
    periods = list(range(1, 4+1))
    
    # Minimum number of fire polygons (imputation threshold)
    imputation_th = 50
    
    # Variables to perform imputation over
    imputation_vars = ['size_q50', 'size_q99', 'season_q50']
    
    
    # Files
    #--------
    
    # Data directory
    in_dir = '../../data/processed/clustering_algorithm/'
    
    # Fire characteristics data file
    data_fp = 'input_vars_grid{}_period{}.csv'
    
    # Edge list file
    edge_fp = 'edge_list_connectivity8_grid_{}_epsg5880.csv'.format(grid)
    
    # Output data file
    out_fp = 'input_vars_grid{}_imputation{}_period{}.csv'
    
    
    
    # -------------------- Data processing ------------------------------------
    
    # Read the edge list
    edges = pd.read_csv(in_dir + edge_fp)
    
    
    # Processing one period at a time
    for p in periods:
        
        # Read data
        data = pd.read_csv(in_dir + data_fp.format(grid, p))
        
        # Get list of cells with <= imputation_th
        imputation_cells = data.loc[data['number_fires'] <= imputation_th, 'id'].tolist()
        
        # For each of these cells
        for cell in imputation_cells:
            
            # Get the neighbours with sufficent fire polygons
            neighbours = get_neighbours(cell, data, edges, imputation_th)
            
            # print(neighbours)
            # print(data.head())
            # print(data.loc[data['id'].isin(neighbours), imputation_vars])
            
            # HINT: there should be an if statement here addressing the possibility 
            #       that there are no enighbours with enough number of fires
            if len(neighbours) == 0:
                print(cell, p)
                
                # In this case, let's try just removing these cells
                data = data.drop(index = data[data['id'] == cell].index)
                
            else:
            
                # Calculate the mean values of imputation_vars among the selected neighbours
                mean_imputation_vars = data.loc[data['id'].isin(neighbours), imputation_vars].mean().tolist()
                
                # Assign the means to the cell
                data.loc[data['id'] == cell, imputation_vars] = mean_imputation_vars
                
                
            
        # Convert season median back to integer (now it is a float)
        data['season_q50'] = data['season_q50'].apply(lambda x: round(x, 0))
            
        # Now that we have updated the values, save the data
        data.to_csv(in_dir + out_fp.format(grid, imputation_th, p), index = False)
    
    
    
    
    
if __name__ == '__main__':
    main()

