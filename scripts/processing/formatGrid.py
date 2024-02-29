# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 11:02:27 2022

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that formats a given grid shapefile by adding the information for:
    - Whether a grid cell is fully within the Cerrado [withinCerr = 0/1]
    - The area of the cell that falls within the Cerrado [clipArea]
    - The percentage of a grid cell within the Cerrado[%area_crop]
    
It requires a grid shapefile as an input where cells have a unique ID. It returns
a grid with columns:
    id | withinCerr | clipArea | %area_crop
    
This script unifies the previous separate scripts that completed the task of 
formatting a grid:
    - labelPolygons_withinPolygon.py
    - calculate_clippedArea.py
    - grid_add_cropArea.py
    
"""

"""
Imports
"""
import geopandas as gpd
import pandas as pd

"""
Functions
"""
def labelCells_withinCerrado(grid, cerrado):
    """
    Function that identifies those grid cells that fall completely within the 
    Cerrado. It creates a new column 'withinCerr' with value 1 if polygons falls 
    completely within the Cerrado, and 0 otherwise.
    
    Performs the same operation as the former script labelPolygons_withinPolygon.py

    Parameters
    ----------
    grid : GeoDataframe
        The grid.
    cerrado : GeoDataframe
        The Cerrado polygon.

    Returns
    -------
    grid : GeoDataframe
        The grid.

    """
    
    # Preparing input ---------------------------------------------------------
    
    # We are going to work on the Cerrado's gdf CRS
    # CRS in which we want to work - we use this as it is the CRS of MapBiomas
    # data
    crs = "EPSG:4326"
    
    # Recording original CRS of grid, to later convert it back
    grid_crs = grid.crs
    cerrado_crs = cerrado.crs
    
    # Converting data to working crs
    grid = grid.to_crs(crs)
    cerrado = cerrado.to_crs(crs)
    
    # Identifying polygons fully within Cerrado -------------------------------
    # Creating column of 0s
    grid['withinCerr'] = 0
    
    # Labelling all of those polygons within with a 1
    grid.loc[grid['geometry'].within(cerrado.loc[0, 'geometry']), 'withinCerr'] = 1

    # Converting CRS back to original
    grid = grid.to_crs(grid_crs)
    cerrado = cerrado.to_crs(cerrado_crs)
    
    
    return grid
    
def calculate_clippedArea(grid, cerrado):
    """
    Function that calculates the area of each cell that is within the Cerrado.
    It does so by clipping each with the Cerrado polygon, and then calculating
    the area of the clipped geometry for each cell.
    
    It casts the information into a new column 'clipArea'.
    
    Performs the same operation as the former script calculate_clippedArea.py
    

    Parameters
    ----------
    grid : GeoDataframe
        The grid.
    cerrado : GeoDataframe
        The Cerrado polygon.

    Returns
    -------
    grid : GeoDataframe
        The grid.

    """
    
    # Preparing input ---------------------------------------------------------
    
    # We are going to work on a projected CRS
    # CRS in which we want to work - we use this as it is the CRS of MapBiomas
    # data
    crs = "EPSG:5880"
    
    # Recording original CRS of grid, to later convert it back
    grid_crs = grid.crs
    cerrado_crs = cerrado.crs
    
    # Converting data to working crs
    grid = grid.to_crs(crs)
    cerrado = cerrado.to_crs(crs)
    
    # Calculating cells' areas within Cerrado ---------------------------------
    
    # First, subset into different GeoDataframe those cells within the Cerrado
    # and the rest
    grid_within = grid[grid['withinCerr'] == 1].reset_index(drop = True)
    grid_not_within = grid[grid['withinCerr'] == 0].reset_index(drop = True)
    
    print('   1. Cells within the Cerrado')
    
    # Then, for cells within the Cerrado, clipArea is just the cell's area
    grid_within['clipArea'] = grid_within['geometry'].apply(lambda x: x.area)
    
    
    # For cells not fully within the Cerrado, we clip the geometries and then 
    # calculate their areas

    print('   2. Cells not fully within the Cerrado')

    # Obtaining the list of clipped cells
    clipped_grid = gpd.clip(grid_not_within['geometry'], cerrado.loc[0, 'geometry'])

    # Calculating areas of clipped cells
    grid_not_within['clipArea'] = clipped_grid.geometry.area
    
    
    # Finally, concatenating both subsets
    grid = gpd.GeoDataFrame(pd.concat([grid_within, grid_not_within], ignore_index = True))
    # and ordering them by index
    grid = grid.sort_values(by = ['id'])
    
    # Converting CRS back to original
    grid = grid.to_crs(grid_crs)
    cerrado = cerrado.to_crs(cerrado_crs)
    
    return grid


def calculate_percentageArea(grid):
    """
    Function that calculates the percentage area within the Cerrado of each cell.
    It does so by using the clipArea information and calculating the area of 
    each cell.
    
    It casts the information into a new column '%area_crop'.
    
    Performs the same operation as the former script grid_add_cropArea.py
    

    Parameters
    ----------
    grid : GeoDataframe
        The grid.

    Returns
    -------
    grid : GeoDataframe
        The grid.

    """
    
    # Preparing input ---------------------------------------------------------
    
    # We are going to work on a projected CRS
    # CRS in which we want to work - we use this as it is the CRS of MapBiomas
    # data
    crs = "EPSG:5880"
    
    # Recording original CRS of grid, to later convert it back
    grid_crs = grid.crs
    
    # Converting data to working crs
    grid = grid.to_crs(crs)
    
    # Calculate each cell's total area
    grid['cellArea'] = grid['geometry'].apply(lambda x: x.area)
    
    # Percentage area of the cell within Cerrado
    grid['%area_crop'] = 100 * grid['clipArea'] / grid['cellArea']
    
    # Dropping column
    grid = grid.drop(columns = ['cellArea'])
    
    # Converting CRS back to original
    grid = grid.to_crs(grid_crs)
    
    return grid
    
def calculateCentroid_latLon(grid):
    
    # # Storing original CRS
    # original_crs = grid.crs
    
    # # Reprojecting grid's CRS to a geographic one
    # grid = grid.to_crs('epsg:4326')
    
    # Obatining the latitude and longitude coordinates
    grid["lon"] = grid.centroid.x
    grid["lat"] = grid.centroid.y
    
    # # Reporjecting grid to original CRS
    # grid = grid.to_crs(original_crs)
    
    return grid

"""
MAIN
"""
def main():
    
    # ---------------------- User inputs -----------------------------------
    # Grid we want to format
    # grid_fp = '../data/cerrado_grid_08deg/cerrado_grid_02deg_cropped.shp'
    grid_fp = '../../data/shapes/cerrado_grid_50km_cropped.shp'
    
    # Geometry for spatial operation
    cerrado_fp = '../../data/shapes/cerrado_1-250000.shp'
    
    # --------------------- Reading input data -----------------------------
    
    grid = gpd.read_file(grid_fp)
    cerrado = gpd.read_file(cerrado_fp)
    
    # --------------------- Formatting grid --------------------------------
    
    # 1. Labelling cells fully within the Cerrado
    print('Labelling polygons within Cerrado')
    grid = labelCells_withinCerrado(grid, cerrado)
    
    # 2. Calculating the area within the Cerrado of each cell
    print('Calculating clip area')
    grid = calculate_clippedArea(grid, cerrado)
    
    # 3. Calculating percentage area of cells within the Cerrado
    print('Calculating percentage area within the Cerrado')
    grid = calculate_percentageArea(grid)
    
    # 4. Calculating the latitude and longitude coordinates in a geographic projection
    print('Calculating lat lon coordinates in a geographic CRS')
    grid = calculateCentroid_latLon(grid)
    
    # 5. Selecting relevant columns
    grid = grid[['id', 'withinCerr', 'clipArea', '%area_crop', 'lon', 'lat', 'geometry']]
    
    
    # --------------------- Overwrite gdf ----------------------------------
    grid.to_file(grid_fp)


if __name__ == "__main__":
    main()

