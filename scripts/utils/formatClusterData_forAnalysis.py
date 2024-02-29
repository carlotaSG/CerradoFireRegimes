# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 17:40:23 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script containing functions that prepare the data frames for various variables
per cluster for analysis.

"""

"""
Imports
"""
import pandas as pd
import geopandas as gpd

# Importing cluster colours and alphabetical codes depending on the method chosen
import sys
sys.path.append('./utils/')
from cluster_identifier import four_periods



"""
Functions
"""

def periods_list(n_periods):
    
    if n_periods == 4:
        periods = ['(1985, 1993)', '(1994, 2002)', '(2003, 2011)', '(2012, 2020)']
        
    return periods

def formatData_ksolutions(nyears, nperiods, pid, k):
    
    # Getting the dictionary relating cells to their cluster, alphabetical code and colour code
    if nperiods == 4:
        clust_ids = four_periods()
    
    # Files with the classification solutions for the period of interest
    in_dir = '../../analysis/{}_year_periods/'.format(nyears)
    
    if pid == 1:
        data = pd.read_csv(in_dir + 'constClusteringClass_ksolutions_period{}.csv'.format(pid))
        
        # Keeping only the columns we are intersted with
        data = data[['id', 'clust_{}'.format(k)]].rename(columns = {'clust_{}'.format(k): 'period{}'.format(pid)})
    
    else:
        data = pd.read_csv(in_dir + 'constClusteringEvol_from{}solutions_afterMerge2_period{}.csv'.format(k, pid))
        
        # Keeping only the columns we are intersted with
        data = data[['id', 'merge_clust']].rename(columns = {'merge_clust': 'period{}'.format(pid)})
    
    # Creating columns with the cluster's letter identifier and corresponding colour
    # Fetching the dictionary of alphabetical-colur codes per cluster
    period_ids = clust_ids['period{}'.format(pid)]
    
    # Creating two columns for this period, one will store the letter idenitfier and the other one the colour code
    data['clustid{}'.format(pid)] = '99'
    data['col{}'.format(pid)] = '#FFFFFF'
    
    # Working on a dictionary element at a time (a cluster)
    for cell in period_ids.keys():
    
        # Fetching the cluster to which this cell belongs in period p
        clust = data.loc[data['id'] == cell, 'period{}'.format(pid)].item()
    
        # All the cells belonging to the same cluster in the period are assigned the same letter
        data.loc[data['period{}'.format(pid)] == clust, 'clustid{}'.format(pid)] = period_ids[cell]['clust_id']
        data.loc[data['period{}'.format(pid)] == clust, 'col{}'.format(pid)] = period_ids[cell]['clust_col']
        
    
    # Adding the grid geometrical information
    grid = gpd.read_file('../../data/shapes/cerrado_grid_50km_cropped.shp')[['id', '%area_crop', 'geometry']]
    
    data = grid.merge(data, on = 'id', how = 'right')
    
    return data


def formatData_expVariables(data, nperiods, nyears, p, clim = True, lu = True, popdens = True,
                            livestock = True, topography = True, fragmentation = True,
                            cell_fragmentation = True):
    
    periods = periods_list(nperiods)
    
    if clim == True:
        # Reading climatic data
        clim_dat = pd.read_csv('../../data/processed/climate_data/summary_period_allVariables_{}periods.csv'.format(nperiods))
        
        # Keeping only period of interest
        clim_dat = clim_dat[clim_dat['period'] == periods[p-1]].reset_index(drop = True)
        clim_dat = clim_dat.drop(columns = 'period')
    
        # Aggregating to clustering classification (along with the geometric data)
        data = data.merge(clim_dat, on = 'id', how = 'left')
        
        del clim_dat
        
    if lu == True:
        # Reading hte land-use and land-use change data
        lu_dat = pd.read_csv('../../data/processed/summary_tables/{}_year_periods/lulc_percArea_allPeriods.csv'.format(nyears))
        
        # Keeping only period of interest
        lu_dat = lu_dat[lu_dat['period'] == periods[p-1]].reset_index(drop = True)
        lu_dat = lu_dat.drop(columns = 'period')
        
        # Aggregating geometric information
        data = data.merge(lu_dat, on = 'id', how = 'left')
        
        del lu_dat
        
    if popdens == True:
        # Reading the population density data
        pop_dat = pd.read_csv('../../data/processed/summary_tables/{}_year_periods/population_density_allPeriods.csv'.format(nyears))
        
        # Keeping only period of interest
        pop_dat = pop_dat[pop_dat['period'] == periods[p-1]].reset_index(drop = True)
        pop_dat = pop_dat.drop(columns = 'period')
        
        # Aggregating geometric information
        data = data.merge(pop_dat, on = 'id', how = 'left')
        
        del pop_dat
        
    if livestock == True:
        # Reading the livestock density data
        lstock_dat = pd.read_csv('../../data/processed/summary_tables/{}_year_periods/livestock_density_allPeriods.csv'.format(nyears))
        # Keeping only period of interest
        lstock_dat = lstock_dat[lstock_dat['period'] == periods[p-1]].reset_index(drop = True)
        lstock_dat = lstock_dat.drop(columns = 'period')
        
        # Aggregating geometric information
        data = data.merge(lstock_dat, on = 'id', how = 'left')
        
        del lstock_dat
        
    if topography == True:
        # reading the topography data
        top_dat = pd.read_csv('../../data/processed/summary_tables/topography_allVars_grid50km.csv')
        
        # Aggregating geometric information
        data = data.merge(top_dat, on = 'id', how = 'left')
        
        del top_dat
        
    if fragmentation == True:
        # reading hte landscape fragmentation measures
        frag_dat = pd.read_csv('../../data/processed/summary_tables/{}_year_periods/landscape_metrics_allPeriods.csv'.format(nyears))
        
        # Aggregating geometric information
        data = data.merge(frag_dat, on = 'id', how = 'left')
        
        del frag_dat
        
    if cell_fragmentation == True:
        # reading hte landscape fragmentation measures
        frag_dat = pd.read_csv('../../data/processed/summary_tables/{}_year_periods/landscape_cell_metrics_allPeriods.csv'.format(nyears))
        
        # Aggregating geometric information
        data = data.merge(frag_dat, on = 'id', how = 'left')
        
        del frag_dat
        
        
    return data
        
    
def unitsLabelsNames_expVariables():
    
    the_info = {
        # Climate
        'clim' : {
            'vars' : ['mean_prec', 'mean_temp', 'mean_rh', 'mean_vpd', 'mean_MCWD_aet', 'mean_prec_ds', 'mean_temp_ds', 'mean_peak_ds'],
            'labels' : ['Mean annual total precipitation', 'Mean annual temperature', 
                        'Mean annual relative humidity', 'Mean annual VPD', 'Mean annual MCWD',
                        'Mean annual total\ndry-season precipitation', 'Mean dry-season temperature',
                        'Peak of the dry season'],
            'units' : ['mm', 'C', '%', 'hPa', 'mm', 'mm', 'C', 'month']
            },
        # Land-use change
        'luch' : {
            'vars' : ['natural_delta', '4_delta', '3_delta', '12_delta', '11_delta', '15_delta', '18_delta'],
            'labels' : ['Change in Natural area', 'Change in Savanna area', 'Change in Forest area', 
                        'Change in Grassland area','Change in Wetland area', 'Change in Pasture area', 
                        'Change in Cropland  area']
            },
        # Land cover composition
        'lu' : {
            'vars' : ['natural', 'anthropic', '4', '3', '12', '11', '15', '18'],
            'labels' : ['Natural area', 'Anthropic area', 'Savanna area', 'Forest area', 
                        'Grassland area','Wetland area', 'Pasture area', 'Cropland  area']
            },
        
        # Population density
        'popdens' : {
            'vars' : ['areaw_popdensity', 'areaw_popdensity_delta'],
            'labels' : ['Population density', 'Change in population density'],
            'units' : ['person/km2', 'person/km2']
            },
        # Livestock density
        'livestock' : {
            'vars' : ['areaw_dens', 'areaw_dens_delta'],
            'labels' : ['Livestock density', 'Change in livestock density'],
            'units' : ['animal/km2', 'animal/km2']
            },
        # Topograhy
        'topography': {
            'vars' : ['roughness', 'slope', 'TRI', 'elevation'],
            'labels' : ['Roughness', 'Slope', 'Topographic Roughness Index', 'Elevation'],
            'units' : ['?', '%', '?', 'm']
            },
        # Landscape fragmentation
        'fragmentation?cell': {
            'vars' : ['proportion_of_landscape', 'number_of_patches', 'patch_density', 
                      'largest_patch_index', 'total_edge', 'edge_density', 
                      'landscape_shape_index', 'area_mn', 'area_md', 'area_cv'],
            'labels': ['Proportion of landscape','Number of patches', 'Patch density', 
                       'Largest patch index', 'Total edges', 'Edge density',
                       'Landscape shape index', 'Mean patch area', 'Median patch area', 'Coeff. of Variation of patch area'],
            'units': ['%', 'n', 'km-2', '?', 'n', 'km-2', '?', 'km2', 'km2', '?']
            },
        # Landscape cell fragmentation
        'fragmentation': {
            'vars' : ['number_of_patches', 'patch_density', 'total_edge', 'edge_density', 
                      'landscape_shape_index', 'area_mn', 'area_cv'],
            'labels': ['Number of patches', 'Patch density', 'Total edges', 'Edge density',
                       'Landscape shape index', 'Mean patch area', 'Coeff. of Variation of patch area'],
            'units': ['n', 'km-2', 'n', 'km-2', '?', 'km2', '?']
            }
        }
    
    return the_info


# FUTNION FOR FIRE CHARACTERISTICS AND THE DICITONARU