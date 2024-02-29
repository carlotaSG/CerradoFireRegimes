# ==============================================================================
#
# Created: 11/01/2022
# Author: Carlota Segura-Garcia
# Email: carlota.seguragarcia@ouce.ox.ac.uk
#
# UPDATED 22/08/2023: to adapt it to work for chapter 02, just changing the file
#                     directories.
#
# Script that calculates different climatic magnitudes from Terraclimate data for each 
# cell in a grid laid over the Cerrado.
#
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

precipitation <- function(grid_month, terrac_fp, year, month){
  # Function that, for a given grid, month, year and data file, reads the precipitation
  # data and calculates certain precipitation-related climate variables for each cell
  
  # Terraclimate Precipitation variable name
  variable <- 'pr'
  
  # Converting month to a string of two digits
  month_str <- sprintf("%02d", month)
  # Data file with precipitation data for year and month
  prec_fp <- glue(terrac_fp)
  #print(prec_fp)
  # Reading tif file
  prec_data <- raster(prec_fp)
  
  # print(crs(prec_data))
  # print(crs(grid_month))

  # # Calculate mean monthly precipitation (units of m/day) for each cell
  a_list <- foreach(i=1:nrow(grid_month), .packages='exactextractr') %dopar% {
    exact_extract(prec_data[[1]], grid_month[i,], 'mean')
    # what is this tpe of data is it a list or a number
  }

  grid_month$mean_prec <- unlist(a_list)
  
  return(grid_month)
}

aet <- function(grid_month, terrac_fp, year, month){
  
  # Terraclimate Actual Evapotranspiration variable name
  variable <- 'aet'
  
  # Converting month to a string of two digits
  month_str <- sprintf("%02d", month)
  # Data file with precipitation data for year and month
  aet_fp <- glue(terrac_fp)
  # print(temp_fp)

  # Reading tif file
  aet_data <- raster(aet_fp)
  
  # Calculate mean monthly precipitation (units of m/day) for each cell
  a_list <- foreach(i=1:nrow(grid_month), .packages='exactextractr') %dopar% {
    exact_extract(aet_data[[1]], grid_month[i,], 'mean')
    # what is this tpe of data is it a list or a number 
  }
  
  grid_month$mean_aet <- unlist(a_list)
  
  # Finally, Terraclimate AET requires a rescaling of 0.1
  grid_month$mean_aet <- grid_month$mean_aet * 0.1
  
  return(grid_month)
}

pet <- function(grid_month, terrac_fp, year, month){
  
  # Terraclimate Potential Evapotranspiration variable name
  variable <- 'pet'
  
  # Converting month to a string of two digits
  month_str <- sprintf("%02d", month)
  # Data file with precipitation data for year and month
  pet_fp <- glue(terrac_fp)
  # print(temp_fp)
  
  # Reading tif file
  pet_data <- raster(pet_fp)
  
  # Calculate mean monthly precipitation (units of m/day) for each cell
  a_list <- foreach(i=1:nrow(grid_month), .packages='exactextractr') %dopar% {
    exact_extract(pet_data[[1]], grid_month[i,], 'mean')
    # what is this tpe of data is it a list or a number 
  }
  
  grid_month$mean_pet <- unlist(a_list)
  
  # Finally, Terraclimate PET requires a rescaling of 0.1
  grid_month$mean_pet <- grid_month$mean_pet * 0.1
  
  return(grid_month)
}

# ==============================================================================
# MAIN

# ----------------- Working directory ------------------------------------------
# Setting working directory to script folder
# Comment out these lines if running on Linux cluster
# current_path = rstudioapi::getActiveDocumentContext()$path 
# setwd(dirname(current_path ))


# ----------------- User inputs ------------------------------------------------

# Grid we are working to retrieve the climatic variables for
grid_size <- '30'
grid_units <- 'km'

# Year period over which to work
start_year <- 1985
end_year <- 2022

# Filepath to datafiles with basename
# terrac_fp <- 'D:/carlota/projects/data/Terraclimate_cerrado/{year}{month_str}.{variable}_cerrado.tif'
terrac_fp <- '../../data/raw/terraclimate/{year}{month_str}.{variable}.tif'

# Output filepath
out_fp <- '../../data/processed/climate_data/terraclimate_grid{grid_size}{grid_units}_{year}.csv'


# ---------------- Reading grid ------------------------------------------------

# Filepath
grid_fp <- glue('../../data/shapes/cerrado_grid_{grid_size}{grid_units}_cropped.shp')
                
# Read the grid into memory
grid <- st_read(grid_fp)
# print(crs(grid))
# print(class(grid))
# Reading one raster file to get CRS
raster_file <- raster('../../data/raw/terraclimate/199001.pr.tif')
# print(crs(raster_file))
grid <- st_transform(grid, crs(raster_file))
# print(crs(grid))

# TODO: check grid columns here and discard useless

# ----------------- Setup cluster ----------------------------------------------
cl <- makeCluster(28)
registerDoParallel(cl)

# ------------------------------------------------------------------------------
# ----------------- Calculating summary variables per year ---------------------

# Working over years consecutively
for (year in seq(start_year, end_year, 1)){

  print(paste('Working on year', year))

  # Creating list to store the dataframes, one for each month, with all the
  # climatic variables
  grid_year <- c()

  # For each month, calculate the different climatic variables
  for (month in seq(1, 12, 1)){
    print(paste('Month', month))
    tic()
    # Create dataframe for this month as a copy of grid - with only relevant columns
    grid_month <- grid %>% select("id", "geometry")
    
    # Add the year as a column
    grid_month['year'] <- year

    # Add the month as a column
    grid_month['month'] <- month


    # 1. Precipitation
    grid_month <- precipitation(grid_month, terrac_fp, year, month)
    
    # 2. Actual Evapotranspiration
    grid_month <- aet(grid_month, terrac_fp, year, month)
    
    # 3. Potential Evapotranspiration
    grid_month <- pet(grid_month, terrac_fp, year, month)

    # Converting monthly information to dataframe
    grid_month <- grid_month %>% st_drop_geometry()

    # Appending grid for this month to list of grids for this year
    grid_year <- append(grid_year, list(grid_month), 0)
    # print(head(grid_month))
    toc()
  }

  # Finally, we stack all dataframes (one per month) and the list into a unique dataframe
  grid_year <- bind_rows(grid_year)
  
  # Temporal modification!!! (05/04/2022)
  # Read the data that we already have for this year
  # data <- read.csv(glue(out_fp))
  
  # data <- merge(x = data, y = grid_year, by = c('id', 'year', 'month'), all = TRUE)

  # Now we can write the object to output file
  write.csv(grid_year, glue(out_fp), row.names = FALSE)

}



# set grud size and range of years to work on

# work in parallel over the different years in the series.

# for each year, loop over the 12 months.

# for each month, calculate the different variables for each cell in grid, one at a time and per theme
# so that tif files are opened and read only once