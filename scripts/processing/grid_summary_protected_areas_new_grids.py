# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 13:59:50 2022

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculate the area fraction per cell in a grid occupied by protected
areas. It does so per year, considering only those protected areas created before
or on the same year. Returns a dataframe per year with the total area fraction
occupied by protected areas, and then by category:
    cell_ID | year | total | UC | TI | PI | US | cat_1 | cat_2 | ...
    
where:
    total = UC + TI
    UC = PI + US        (Unidade de Conservacao)
    PI                  (Protecao inegral)
    US                  (Uso Sustentavel)
    TI                  (Terras indigenas)
    cat_i               (Type of Protected Area)
    

"""


"""
Imports
"""
import geopandas as gpd
import pandas as pd


import sys
sys.path.append('../../scripts')
import fileList
import calculate_fracAreaClipped as fracArea

"""
Functions
"""
def group_and_sum(data, group):
    # Function that groupsby 'id' and group and adds the corresponding rows
    # and calcualtes the total frac_area per group.
    # Then, it formats the dataframe from id | group | frac_area
    # to id | group1 | group2
    
    # Groupby and summarise
    data_group = data[['id', group, 'frac_area']].groupby(by = ['id', group]).sum().reset_index()
    
    # Reformat from long to wide format
    data_group = pd.pivot(data_group, index = 'id', columns = group, values = 'frac_area').reset_index()
    
    # Filling NAs with 0s
    data_group = data_group.fillna(0)
    
    return data_group
    



"""
MAIN
"""
if __name__ == '__main__':
    
    # -------------------------------------------------------------------------
    # User inputs
    
    # Shapefile with Protected areas overlapping the Cerrado
    pa_fp = '../../data/raw/protected_areas/protected_areas_cerrado_fixed.shp'
    
    # Grid SHP
    grid_id = '30km'
    grid_fp = '../../data/shapes/cerrado_grid_{}_cropped.shp'.format(grid_id)
    
    # Threshold to select cells
    area_thresh = 0.0 #75.0
    
    # include_cells = [567, 568]
    # exclude_cells = [1451]

    
    # Years to work with
    start_year, end_year = 1985, 2020
    
    # Output filepath
    out_fp = '../../data_new_grids/processed/protected_areas/protected_areas_{}_{}.csv'
    
    
    # -------------------------------------------------------------------------
    # Calculating the percentage area
    
    # First, let's read the grid and the protected area's information
    # Reading the grid and subsetting relevant cells
    grid = gpd.read_file(grid_fp)
    # grid = grid[(grid['%area_crop'] >= area_thresh) | grid['id'].isin(include_cells)]
    # grid = grid[~grid['id'].isin(exclude_cells)].reset_index(drop = True)
    grid = grid.to_crs('epsg:5880')
    
    pa = gpd.read_file(pa_fp)
    pa = pa.to_crs('epsg:5880')
    
    

    
    # Working on the years consecutively
    for year in range(start_year, end_year + 1):
        
        print('Working on year {}'.format(year))
        
        # Select those protected areas created on or before year
        pa_year = pa[pa['cre_year'] <= year]
        
        # print(pa_year)
        
        part_1 = pa_year[['category', 'geometry']].dissolve(by = 'category').reset_index()
        part_2 = pa_year[['pa_type', 'geometry']].dissolve(by = 'pa_type').reset_index()
        part_2 = part_2.rename(columns = {'pa_type': 'category'})
        
        # For the total area
        pa_year['pa_total'] = 'pa_total'
        part_3 = pa_year[['pa_total', 'geometry']].dissolve(by = 'pa_total').reset_index()
        part_3 = part_3.rename(columns = {'pa_total': 'category'})
        
        # For the Conservation Units (not TI)
        part_4 = pa_year.loc[pa_year['pa_type'] != 'TI', ['pa_type', 'geometry']]
        part_4['category'] = 'CU'
        part_4 = part_4.dissolve(by = 'category').reset_index()
        # part_4 = part_4.rename(columns = {'pa_type': 'category'})
        # print(part_4)
        
        data_year = gpd.GeoDataFrame(pd.concat([part_1, part_2, part_3, part_4], ignore_index = True), crs = part_1.crs)
        # Erasing duplicate TI
        data_year = data_year.drop_duplicates(subset = 'category')
        
        # print(data_year)
        
        # Freeing memory
        del pa_year, part_1, part_2, part_3, part_4
        
        # Calculate cell area fraction occupied per protected area category
        frac_areas = fracArea.clip_calculateArea_withPool(data_year, grid, 'category', 'id', 10)
        
        # Rename geom column to 'id'
        frac_areas = frac_areas.rename(columns = {'geom':'id'})
        
        # print(frac_areas)
        
        frac_areas = pd.pivot(frac_areas, index = 'id', columns = 'category', values = 'frac_area').reset_index()
        
        # Fill NAs with 0s
        frac_areas = frac_areas.fillna(0)
        
        # Insert year column
        frac_areas.insert(1, 'year', year)
        
        # print(frac_areas)
        # break
    
        # Saving file
        print('Saving protected areas file for year {}.\n'.format(year))
        frac_areas.to_csv(out_fp.format(grid_id, year), index = False)
        
        # Freeing memory
        del frac_areas
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        # # Calculate cell area fraction occupied per protected area
        # frac_areas = fracArea.clip_calculateArea_withPool(pa_year, grid, 'pa_id', 'id', 10)
        
        # # Rename geom column to 'id'
        # frac_areas = frac_areas.rename(columns = {'geom':'id'})
        
        # # Merge the information of each protected area
        # frac_areas = frac_areas.merge(pa_year, on = 'pa_id', how = 'left')
        
        # # Aggregate by pa_type and format
        # frac_areas_pa = group_and_sum(frac_areas, 'pa_type')
        
        # print(frac_areas_pa)
        
        # # Aggregate by category and format
        # frac_areas_cat = group_and_sum(frac_areas[frac_areas['pa_type'] != 'TI'].reset_index(), 'category')
        
        # print(frac_areas_cat)
        
        # # Merge back to grid (converted to geometry): reason is that not all cells 
        # # overlap with protected areas, these should be assigned a value 0
        # grid_frac_area = grid[['id']].merge(frac_areas_pa, on = ['id'], how = 'left')
        # grid_frac_area = grid_frac_area.merge(frac_areas_cat, on = ['id'], how = 'left')
        
        # # Filling NAs with 0s
        # grid_frac_area = grid_frac_area.fillna(0)
        
        # TODO: calculate overlapping geometries in the year and calculate the fraction of area 
        # of each cell that they occupy
        
        # # Alternative manner of finding overlapping polygons
        # sjoin_gdf = gpd.sjoin(pa_year[['pa_']], pa_year)
        # print(sjoin_gdf[['pa_id_right', 'pa_id_left', 'index_right']])
        
    
        
        # sjoin_gdf = sjoin_gdf.loc[sjoin_gdf.index != sjoin_gdf.index_right]
        
        # print(sjoin_gdf[['pa_id_right', 'pa_id_left']])
        
        # # creating column with set(pa_id_right, pa_id_left)
        # sjoin_gdf['sets'] = list(zip(sjoin_gdf['pa_id_right'], sjoin_gdf['pa_id_left']))
        # sjoin_gdf['sets'] = sjoin_gdf['sets'].apply(lambda x: frozenset(x))
        
        # print(sjoin_gdf[['pa_id_right', 'pa_id_left', 'sets']])
        
        # sjoin_gdf = sjoin_gdf.drop_duplicates(subset = 'sets', ignore_index = True)
        
        # print(sjoin_gdf[['pa_id_right', 'pa_id_left', 'sets']])
        # print(sjoin_gdf.columns)
        
        # If I want to have area fraction per subcategory of protected area,
        # then it is quite complicated because protected areas of various categories
        # Are overlapping with one another... and so discounting the overlap is
        # complicated to achieve. 
        
        # I could calculate the area fraction in a different strategy:
        # For the subcategories:
            # Dissolve all polygons converting into a single Multipolygon,
            # calculate the total area fraction occupied by this multipolygon in each cell
            # by those polygons
        # DO the same for the categories and for the total
        
        
        # # Summarise:
        # # TODO: insert these columns in a sepcific position
        # grid_frac_area['CU'] = grid_frac_area['PI'] + grid_frac_area['US']
        # grid_frac_area['total'] = grid_frac_area['CU'] + grid_frac_area['TI']
        
        # # Add year column
        # grid_frac_area['year'] = year
        
        # print(grid_frac_area)
        
        # break
            

        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
