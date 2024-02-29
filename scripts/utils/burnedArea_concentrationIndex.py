# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 13:42:31 2022

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculates Burned area precipitation index as in paper Defining pyromes, inspired
on Seasonality of Precipitation in the US.
"""

"""
Imports
"""
import pandas as pd
import numpy as np

"""
Functions
"""
def month_radian_equivalence():
    
    equivalence = pd.DataFrame(
        data = {
            'month' : range(1, 12 + 1),
            'rad'   : np.arange(15, 345 + 30, 30) * np.pi / 180
            }
        )
    
    return equivalence
    
def month_to_radians(months):
    # Gets a list of months as a Pandas Series and returns it as a list of angles
    # expressed in radians
    
    # monthsToRadians = {
    #     '1' : 15 * np.pi/180,
    #     '2' : 45 * np.pi/180,
    #     '3' : 75 * np.pi/180,
    #     '4' : 105 * np.pi/180,
    #     '5' : 135 * np.pi/180,
    #     '6' : 165 * np.pi/180,
    #     '7' : 195 * np.pi/180,
    #     '8' : 225 * np.pi/180,
    #     '9' : 255 * np.pi/180,
    #     '10': 285 * np.pi/180,
    #     '11': 315 * np.pi/180,
    #     '12': 345 * np.pi/180
    #     }
    
    eq = month_radian_equivalence()
    # print(eq)
    # print(months)
    # Converting every month to a radian
    # months_rad = [eq.loc[eq['month'] == m, 'rad'] for m in months]
    
    months_rad = months.apply(lambda m: eq.loc[eq['month'] == m, 'rad'].item())
    # print(months_rad)
    
    return months_rad
    
def polar_to_cartesian(ba_polar):
    # Function transforming polar coordantes (angle, magnitude):(month_rad, ba)
    # to cartesian coordiantes (x, y)
    
    # X coordinate
    x = ba_polar['ba'] * np.cos(ba_polar['month_rad'])
    # Y coordinate
    y = ba_polar['ba'] * np.sin(ba_polar['month_rad'])
    
    return x, y

def cartesian_to_radians(x, y):
    # Function that converts a point in cartesian to polar coordinates with angle
    # in radians in the range [0, 2*np.pi]
    
    magnitude = np.sqrt(x**2 + y**2)
    # Angle in the range [-np.pi, np.pi]
    angle_deg = np.arctan2(y, x)
    
    # Convert angle to (0, 360)
    if angle_deg < 0:
        angle_deg = angle_deg + 2 * np.pi
    
    return angle_deg, magnitude


def rad_to_month(val_rad):
    # Function that converts a degree to the closest month.
    
    # First, declaring the Data Frame with the degree - month equivalence
    eq = month_radian_equivalence()
    
    # Then, find the month corresponding the the angle closest to val_rad
    val_rad_exact = eq.loc[eq['rad'].sub(val_rad).abs().idxmin(), 'month']
    
    return val_rad_exact

def burnedArea_concentrationIndex(ba_month, ba_col):
    # main function recevieng data frame of month | year | burned_area for a 
    # certain year
    # print(ba_month)
    # 0. Rename the column containing burned area to 'ba'
    ba_month = ba_month.rename(columns = {ba_col : 'ba'})
    # print(ba_month)
    # 1. Converting months to angles (in radians)
    ba_month['month_rad'] = month_to_radians(ba_month['month'])
    
    # print(ba_month)
    
    # 2. Convert pairs month_rad, burned area (angle, magnitude) in polar coordinates
    #    to cartesian coordinates
    ba_month['x'], ba_month['y'] = polar_to_cartesian(ba_month[['month_rad', 'ba']])
    
    
    # 3. Adding vectors in cartesian coordinates
    #    Resulting in a point in cartesian coordinates
    x_period = ba_month['x'].sum()
    y_period = ba_month['y'].sum()
    
    # 4. Convert from cartesian to radians
    month_rad, ba_magnitude = cartesian_to_radians(x_period, y_period)
    
    # 5. Convert radians to month to have the peak of the season
    peak_month = rad_to_month(month_rad)
    
    # 6. Finally, calculate the Burned Area Concentration Index
    ba_ci = 100 * ba_magnitude / ba_month['ba'].sum()
    
    return ba_ci, peak_month
    
    
    
    
    
    
    