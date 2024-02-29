# ==============================================================================
#
# Created: 19/04/2023
# Author: Carlota Segura-Garcia
# Email: carlota.seguragarcia@ouce.ox.ac.uk
#
# Script that calculates average monthly precipitation using CHIRPS data for each 
# cell in a grid laid over the Cerrado.
#
# Sligthly modified from the script ./_chapter01/scripts/grid_summary_chirps.R
#
# ==============================================================================

# ==============================================================================
# Imports
library(raster)
library(sf)
library(exactextractr)
library(rgeos)
library(Rcpp)
library(rgdal)
library(dplyr)
library(tidyverse)
library(doParallel)
library(glue)
library(tictoc)


# ==============================================================================
# Functions

precipitation <- function(grid_month, chirps_fp, year, month){
  # Function that, for a given grid, month, year and data file, reads the precipitation
  # data and calculates certain precipitation-related climate variables for each cell
  
  # CHIRPS Precipitation variable name
  variable <- 'precipitation'
  
  # Converting month to a string of two digits
  month_str <- sprintf("%02d", month)
  # Data file with precipitation data for year and month
  prec_fp <- glue(chirps_fp)

  # Reading tif file
  prec_data <- raster(prec_fp)
  
  # Calculate mean monthly precipitation (units of m/day) for each cell
  a_list <- foreach(i=1:nrow(grid_month), .packages='exactextractr') %dopar% {
    exact_extract(prec_data[[1]], grid_month[i,], 'mean')
    # what is this tpe of data is it a list or a number 
  }

  grid_month$mean_prec <- unlist(a_list)
  
  # CHIRPS units are already mm/month
  
  return(grid_month)
}



# ==============================================================================
# MAIN

# ----------------- Working directory ------------------------------------------
# Setting working directory to script folder
# current_path = rstudioapi::getActiveDocumentContext()$path
# setwd(dirname(current_path ))

print(getwd())
# ----------------- User inputs ------------------------------------------------

# Grid we are working to retrieve the climatic variables for
grid_size <- '30'
grid_units <- 'km'

# Year period over which to work
start_year <- 1985
end_year <- 2022

# Filepath to datafiles with basename
# chirps_fp <- 'D:/carlota/projects/data/CHIRPS_cerrado/{year}{month_str}.precipitation_cerrado.tif'
chirps_fp <- '../../data/raw/CHIRPS_cerrado/{year}{month_str}.precipitation.tif'

# Output filepath
out_fp <- '../../data/processed/climate_data/chirps_grid{grid_size}{grid_units}_{year}.csv'


# ---------------- Reading grid ------------------------------------------------

# Filepath
grid_fp <- glue('../../data/shapes/cerrado_grid_{grid_size}{grid_units}_cropped.shp')
                
# Read the grid into memory
grid <- st_read(grid_fp)

# TODO: check grid columns here and discard useless

# ----------------- Setup cluster ----------------------------------------------
cl <- makeCluster(28)
registerDoParallel(cl)

# ------------------------------------------------------------------------------
# ----------------- Calculating summary variables per year ---------------------

# Working over years consecutively
for (year in seq(start_year, end_year, 1)){
  
  tic()
  print(paste('Working on year', year))
  
  # Creating list to store the dataframes, one for each month, with precipitation
  grid_year <- c()
  
  # For each month, calculate the different climatic variables
  for (month in seq(1, 12, 1)){
    print(paste('Month', month))
    
    # Create dataframe for this month as a copy of grid - with only relevant columns
    grid_month <- grid %>% select("id", "geometry")
    
    # Add the year as a column
    grid_month['year'] <- year
    
    # Add the month as a column
    grid_month['month'] <- month
    
    # 1. Precipitation
    grid_month <- precipitation(grid_month, chirps_fp, year, month)
    
    # Converting monthly information to dataframe
    grid_month <- grid_month %>% st_drop_geometry()
    
    # Appending grid for this month to list of grids for this year
    grid_year <- append(grid_year, list(grid_month), 0)
    # print(head(grid_month))
    
  }
  
  
  # Finally, we stack all dataframes (one per month) and the list into a unique dataframe
  grid_year <- bind_rows(grid_year) 

  # Now we can write the object to output file
  write.csv(grid_year, glue(out_fp), row.names = FALSE)
  
  toc()
  
}