# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 16:18:51 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk


Dictionaries where reference cells that do not change clusters throughout the 
periods get assigned a cluster number/letter and a colour.

"""

"""
Imports
"""
import matplotlib.pyplot as plt
import matplotlib.colors as mcl
import geopandas as gpd

def palette_cluster():
    
    # List the set of colrs in palette
    cmap = plt.get_cmap('tab20')
    list_colors = [mcl.rgb2hex(cmap(i)) for i in range(cmap.N)]
    
    palette = {
        'A' : list_colors[0],
        'B' : list_colors[1],
        'C' : list_colors[2],
        'D' : list_colors[3],
        'E' : list_colors[4],
        'F' : list_colors[5],
        'G' : list_colors[6],
        'H' : list_colors[7],
        'I' : list_colors[8],
        'J' : list_colors[9],
        'K' : list_colors[10],
        'L' : list_colors[11],
        'M' : list_colors[12],
        'N' : list_colors[13]
        }
    
    return palette
    

def four_periods():
    
    colors = palette_cluster()
    
    # Dictionary of the format: cell_id : {cluster_letter, cluster_colour}
    cluster_structure = { 
        'period1' : {
            297  : {'clust_id' : 'A', 'clust_col' : colors['A']},
            734  : {'clust_id' : 'B', 'clust_col' : colors['B']},
            1122 : {'clust_id' : 'C', 'clust_col' : colors['C']},
            1983 : {'clust_id' : 'D', 'clust_col' : colors['D']},
            2108 : {'clust_id' : 'E', 'clust_col' : colors['E']},
            2329 : {'clust_id' : 'F', 'clust_col' : colors['F']},
            2065 : {'clust_id' : 'G', 'clust_col' : colors['G']},                         # Until here these cells serve as reference for all periods
            185  : {'clust_id' : 'H', 'clust_col' : colors['H']},                         # Used in periods 1, 2 and 3
            584  : {'clust_id' : 'I', 'clust_col' : colors['I']},                         # Used in periods 1 and 2
            2020 : {'clust_id' : 'J', 'clust_col' : colors['J']},                         # Used in periods 1, 2, and 4
            1858 : {'clust_id' : 'K', 'clust_col' : colors['K']},
            1575 : {'clust_id' : 'L', 'clust_col' : colors['L']}                          # used in periods 1 and 3
            },
        'period2' : {
            297  : {'clust_id' : 'A', 'clust_col' : colors['A']},
            734  : {'clust_id' : 'B', 'clust_col' : colors['B']},
            1122 : {'clust_id' : 'C', 'clust_col' : colors['C']},
            1983 : {'clust_id' : 'D', 'clust_col' : colors['D']},
            2108 : {'clust_id' : 'E', 'clust_col' : colors['E']},
            2329 : {'clust_id' : 'F', 'clust_col' : colors['F']},
            2065 : {'clust_id' : 'G', 'clust_col' : colors['G']},                         # Until here these cells serve as reference for all periods
            185  : {'clust_id' : 'H', 'clust_col' : colors['H']},                         # Used in periods 1, 2 and 3
            584  : {'clust_id' : 'I', 'clust_col' : colors['I']},                         # Used in periods 1 and 2
            2020 : {'clust_id' : 'J', 'clust_col' : colors['J']},                         # Used in periods 1, 2, and 4
            1155 : {'clust_id' : 'M', 'clust_col' : colors['M']}                          # Used in periods 2, 3, 4
            },
        'period3' : {
            297  : {'clust_id' : 'A', 'clust_col' : colors['A']},
            734  : {'clust_id' : 'B', 'clust_col' : colors['B']},
            1122 : {'clust_id' : 'C', 'clust_col' : colors['C']},
            1983 : {'clust_id' : 'D', 'clust_col' : colors['D']},
            2108 : {'clust_id' : 'E', 'clust_col' : colors['E']},
            2329 : {'clust_id' : 'F', 'clust_col' : colors['F']},
            2065 : {'clust_id' : 'G', 'clust_col' : colors['G']},                         # Until here these cells serve as reference for all periods
            185  : {'clust_id' : 'H', 'clust_col' : colors['H']},                         # Used in periods 1, 2 and 3
            1155 : {'clust_id' : 'M', 'clust_col' : colors['M']},                         # Used in periods 2, 3, 4
            1575 : {'clust_id' : 'L', 'clust_col' : colors['L']}                          # used in periods 1 and 3
            },
        'period4' : {
            297  : {'clust_id' : 'A', 'clust_col' : colors['A']},
            734  : {'clust_id' : 'B', 'clust_col' : colors['B']},
            1122 : {'clust_id' : 'C', 'clust_col' : colors['C']},
            1983 : {'clust_id' : 'D', 'clust_col' : colors['D']},
            2108 : {'clust_id' : 'E', 'clust_col' : colors['E']},
            2329 : {'clust_id' : 'F', 'clust_col' : colors['F']},
            2065 : {'clust_id' : 'G', 'clust_col' : colors['G']},                         # Until here these cells serve as reference for all periods
            2020 : {'clust_id' : 'J', 'clust_col' : colors['J']},                         # Used in periods 1, 2, and 4
            1155 : {'clust_id' : 'M', 'clust_col' : colors['M']}                          # Used in periods 2, 3, 4
            }
        }
    
    return cluster_structure

def two_periods_clusterEvolution():
    
    colors = palette_cluster()
    
    # Dictionary of the format: cell_id : {cluster_letter, cluster_colour}
    cluster_structure = { 
        'period1' : {
            297  : {'clust_id' : 'A', 'clust_col' : colors['A']},
            734  : {'clust_id' : 'B', 'clust_col' : colors['B']},
            1122 : {'clust_id' : 'C', 'clust_col' : colors['C']},
            1983 : {'clust_id' : 'D', 'clust_col' : colors['D']},
            2108 : {'clust_id' : 'E', 'clust_col' : colors['E']},
            2329 : {'clust_id' : 'F', 'clust_col' : colors['F']},
            1525 : {'clust_id' : 'G', 'clust_col' : colors['G']},
            185  : {'clust_id' : 'H', 'clust_col' : colors['H']},
            2020 : {'clust_id' : 'J', 'clust_col' : colors['J']},
            1858 : {'clust_id' : 'K', 'clust_col' : colors['K']},
            1575 : {'clust_id' : 'L', 'clust_col' : colors['L']}
            },
        'period2' : {
            297  : {'clust_id' : 'A', 'clust_col' : colors['A']},
            1122 : {'clust_id' : 'C', 'clust_col' : colors['C']},
            1983 : {'clust_id' : 'D', 'clust_col' : colors['D']},
            2108 : {'clust_id' : 'E', 'clust_col' : colors['E']},
            2329 : {'clust_id' : 'F', 'clust_col' : colors['F']},
            1525 : {'clust_id' : 'G', 'clust_col' : colors['G']},
            2020 : {'clust_id' : 'J', 'clust_col' : colors['J']},
            1155 : {'clust_id' : 'M', 'clust_col' : colors['M']}
            
            }
        }
    
    return cluster_structure


def two_periods_choosingk():
    
    colors = palette_cluster()
    
    # Dictionary of the format: cell_id : {cluster_letter, cluster_colour}
    cluster_structure = { 
        'period1' : {
            297  : {'clust_id' : 'A', 'clust_col' : colors['A']},
            734  : {'clust_id' : 'B', 'clust_col' : colors['B']},
            1122 : {'clust_id' : 'C', 'clust_col' : colors['C']},
            1983 : {'clust_id' : 'D', 'clust_col' : colors['D']},
            2108 : {'clust_id' : 'E', 'clust_col' : colors['E']},
            2329 : {'clust_id' : 'F', 'clust_col' : colors['F']},
            1525 : {'clust_id' : 'G', 'clust_col' : colors['G']},
            185  : {'clust_id' : 'H', 'clust_col' : colors['H']},
            2020 : {'clust_id' : 'J', 'clust_col' : colors['J']},
            1858 : {'clust_id' : 'K', 'clust_col' : colors['K']},
            1575 : {'clust_id' : 'L', 'clust_col' : colors['L']}
            },
        'period2' : {
            297  : {'clust_id' : 'A', 'clust_col' : colors['A']},
            1122 : {'clust_id' : 'C', 'clust_col' : colors['C']},
            1983 : {'clust_id' : 'D', 'clust_col' : colors['D']},
            2329 : {'clust_id' : 'F', 'clust_col' : colors['F']},
            1525 : {'clust_id' : 'G', 'clust_col' : colors['G']},
            1364 : {'clust_id' : 'J', 'clust_col' : colors['J']},
            1858 : {'clust_id' : 'K', 'clust_col' : colors['K']},
            1575 : {'clust_id' : 'L', 'clust_col' : colors['L']},
            1155 : {'clust_id' : 'M', 'clust_col' : colors['M']},
            750  : {'clust_id' : 'N', 'clust_col' : colors['N']}
            }
        }
    
    return cluster_structure


def format_cluster_classification(data, p, nperiods = 4, evol = True, geom_info = True):
    # Function that, given a cell classificatioin with information: id | clust{p}
    # where id is the cell's id, and clust{p} is the cluster numerical classification for period p,
    # it creates two additional columns containing the alhpabetical identifier for
    # the cluster and the corresponding colour. If requested, it also adds a column with the 
    # georeferrenced geometry of each cell.
    
    if nperiods == 4:  
        # Fetching the cell - cluster_id - cluster_colour for this period
        period_ids = four_periods()['period{}'.format(p)]
        
        # Creating two columns for cluster's alphabetical and coloured ids
        data['clustid{}'.format(p)] = '99'
        data['col{}'.format(p)] = '#FFFFFF'
        
        # Casting the corresponding cluster information one cluster at a time
        for cell in period_ids.keys():
            # Fetching the cluster to which this cell belongs in period p
            clust = data.loc[data['id'] == cell, 'period{}'.format(p)].item()
        
            # All the cells belonging to the same cluster in the period are assigned the same letter
            data.loc[data['period{}'.format(p)] == clust, 'clustid{}'.format(p)] = period_ids[cell]['clust_id']
            data.loc[data['period{}'.format(p)] == clust, 'col{}'.format(p)] = period_ids[cell]['clust_col']
            
    
    if geom_info == True:
        # Adding the grid geometrical information
        grid = gpd.read_file('../../data/shapes/cerrado_grid_50km_cropped.shp')[['id', '%area_crop', 'geometry']]
        
        data = grid.merge(data, on = 'id', how = 'right')
        
    return data


def four_periods_boxplotGroups():
    
    groups = [
        ['H', 'A', 'M'], ['B'], ['C','I'], ['K','J','D'], ['L','G','F'], ['E']
        ]
    
    return groups

def two_periods_clusterEvolution_boxplotGroups():
    
    groups = [
        ['H', 'A', 'M', 'B'], ['C'], ['D'], ['J'], ['L','G','K'], ['F'], ['E']
        ]
    
    return groups

def two_periods_choosingK_boxplotGroups():
    
    groups = [
        ['H', 'A', 'M'], ['C','N'], ['B','D'], ['J'], ['L','G'],['K'], ['F','E']
        ]
    
    return groups