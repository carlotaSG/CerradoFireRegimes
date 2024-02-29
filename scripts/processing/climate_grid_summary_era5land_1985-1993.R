# ==============================================================================
#
# Created: 11/01/2022 Author: Carlota Segura-Garcia Email:
# carlota.seguragarcia@ouce.ox.ac.uk
#
# UPDATED 22/08/2023: to adapt it to work for chapter 02, just changing the file
#                     directories.
#
# Script that calculates different climatic magnitudes from ERA 5 data for each
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

precipitation <- function(grid_month, era_fp, year, month){
  # Function that, for a given grid, month, year and data file, reads the precipitation
  # data and calculates certain precipitation-related climate variables for each cell
  
  # ERA 5 Land Precipitation variable name
  variable <- 'total_precipitation'
  
  # Converting month to a string of two digits
  month_str <- sprintf("%02d", month)
  # Data file with precipitation data for year and month
  prec_fp <- glue(era_fp)
  #print(prec_fp)
  # Reading tif file
  prec_data <- raster(prec_fp)
  
  # Calculate mean monthly precipitation (units of m/day) for each cell
  a_list <- foreach(i=1:nrow(grid_month), .packages='exactextractr') %dopar% {
    exact_extract(prec_data[[1]], grid_month[i,], 'mean')
    # what is this tpe of data is it a list or a number 
  }

  grid_month$mean_prec <- unlist(a_list)

  # Some cells do not overlap with the data, fill the corresponding rows with 9999
  # grid_month[is.na(grid_month)] <- 9999
  # print(grid_month)
  # print(sapply(grid_month, class))
  # When converting units to per month instead of day, we take into account that
  # certain months have 30, 31 or 28 days? TODO
  # factor <- 30
  if (month %in% c(1, 3, 5 ,7, 8, 10, 12)){
    factor <- 31
  }else if (month %in% c(4, 6, 9, 11)){
    factor <- 30
  }else if (month == 2){
    factor <- 28
  }
  
  # Finally, convert units from m/day to mm/month
  grid_month$mean_prec <- grid_month$mean_prec * 1000 * factor
  
  return(grid_month)
}

temperature <- function(grid_month, era_fp, year, month){
  
  # ERA 5 Land Temperature variable name
  variable <- 'temperature_2m'
  
  # Converting month to a string of two digits
  month_str <- sprintf("%02d", month)
  # Data file with precipitation data for year and month
  temp_fp <- glue(era_fp)
  # print(temp_fp)

  # Reading tif file
  temp_data <- raster(temp_fp)
  # Setting the nodata value to 9999
  NAvalue(temp_data) <- 9999
  
  # Calculate mean monthly precipitation (units of m/day) for each cell
  a_list <- foreach(i=1:nrow(grid_month), .packages='exactextractr') %dopar% {
    exact_extract(temp_data[[1]], grid_month[i,], 'mean')
    # what is this tpe of data is it a list or a number 
  }
  
  grid_month$mean_temp <- unlist(a_list)
 
  # Finally, convert units from K to Celsius
  grid_month$mean_temp <- grid_month$mean_temp - 273.15
  
  return(grid_month)
}

RH <- function(grid_month, era_fp, year, month){
  
  # ERA 5 Land Relative Humidity variable name
  variable <- 'RH'
  
  # Converting month to a string of two digits
  month_str <- sprintf("%02d", month)
  # Data file with precipitation data for year and month
  rh_fp <- glue(era_fp)
  # print(temp_fp)
  
  # Reading tif file
  rh_data <- raster(rh_fp)
  
  # Calculate mean monthly precipitation (units of m/day) for each cell
  a_list <- foreach(i=1:nrow(grid_month), .packages='exactextractr') %dopar% {
    exact_extract(rh_data[[1]], grid_month[i,], 'mean')
    # what is this tpe of data is it a list or a number 
  }
  
  grid_month$mean_rh <- unlist(a_list)
  
  # Finally, convert units from K to Celsius
  grid_month$mean_rh <- grid_month$mean_rh
  
  return(grid_month)
}

VPD <- function(grid_month, era_fp, year, month){
  
  # ERA 5 Land RElative Humidity variable name
  variable <- 'VPD'
  
  # Converting month to a string of two digits
  month_str <- sprintf("%02d", month)
  # Data file with precipitation data for year and month
  vpd_fp <- glue(era_fp)
  # print(temp_fp)
  
  # Reading tif file
  vpd_data <- raster(vpd_fp)
  
  # Calculate mean monthly precipitation (units of m/day) for each cell
  a_list <- foreach(i=1:nrow(grid_month), .packages='exactextractr') %dopar% {
    exact_extract(vpd_data[[1]], grid_month[i,], 'mean')
    # what is this tpe of data is it a list or a number 
  }
  
  grid_month$mean_vpd <- unlist(a_list)
  
  # Finally, convert units from K to Celsius
  grid_month$mean_vpd <- grid_month$mean_vpd
  
  return(grid_month)
}

potentialET <- function(grid_month, era_fp, year, month){
  
  # ERA 5 Land Potential Evapotranspiration variable name
  variable <- 'potential_evaporation'
  
  # Converting month to a string of two digits
  month_str <- sprintf("%02d", month)
  # Data file with precipitation data for year and month
  pet_fp <- glue(era_fp)
  # print(temp_fp)
  
  # Reading tif file
  pet_data <- raster(pet_fp)
  # Setting the nodata value to 9999
  NAvalue(pet_data) <- 9999
  
  # Calculate mean monthly potential evapotranspiration (units of m/day) for each cell
  a_list <- foreach(i=1:nrow(grid_month), .packages='exactextractr') %dopar% {
    exact_extract(pet_data[[1]], grid_month[i,], 'mean')
    # what is this tpe of data is it a list or a number 
  }
  
  grid_month$mean_pet <- unlist(a_list)
  
  # When converting units to per month instead of day, we take into account that
  # certain months have 30, 31 or 28 days? TODO
  # factor <- 30
  if (month %in% c(1, 3, 5 ,7, 8, 10, 12)){
    factor <- 31
  }else if (month %in% c(4, 6, 9, 11)){
    factor <- 30
  }else if (month == 2){
    factor <- 28
  }
  
  # Finally, convert units from m/day to mm/month
  grid_month$mean_pet <- grid_month$mean_pet * 1000 * factor
  
  # Finally, convert units from K to Celsius
  grid_month$mean_pet <- grid_month$mean_pet
  
  return(grid_month)
}

actualET <- function(grid_month, era_fp, year, month){
  
  # ERA 5 Land Actual Evapotranspiration (AET) variable name
  variable <- 'total_evaporation'
  
  # Converting month to a string of two digits
  month_str <- sprintf("%02d", month)
  # Data file with precipitation data for year and month
  aet_fp <- glue(era_fp)
  # print(temp_fp)
  
  # Reading tif file
  aet_data <- raster(aet_fp)
  # Setting the nodata value to 9999
  NAvalue(aet_data) <- 9999
  
  # Calculate mean monthly actual evapotranspiration (units of m/day) for each cell
  a_list <- foreach(i=1:nrow(grid_month), .packages='exactextractr') %dopar% {
    exact_extract(aet_data[[1]], grid_month[i,], 'mean')
    # what is this tpe of data is it a list or a number 
  }
  
  grid_month$mean_aet <- unlist(a_list)
  
  # When converting units to per month instead of day, we take into account that
  # certain months have 30, 31 or 28 days? TODO
  # factor <- 30
  if (month %in% c(1, 3, 5 ,7, 8, 10, 12)){
    factor <- 31
  }else if (month %in% c(4, 6, 9, 11)){
    factor <- 30
  }else if (month == 2){
    factor <- 28
  }
  
  # Finally, convert units from m/day to mm/month
  grid_month$mean_aet <- grid_month$mean_aet * 1000 * factor
  
  # Finally, convert units from K to Celsius
  grid_month$mean_aet <- grid_month$mean_aet
  
  return(grid_month)
}


# ==============================================================================
# MAIN

# ----------------- Working directory ------------------------------------------
# Setting working directory to script folder
# current_path = rstudioapi::getActiveDocumentContext()$path 
# setwd(getSrcDirectory(function(){})[1])
# setwd(dirname(current_path ))

print(getwd())

# ----------------- User inputs ------------------------------------------------

# Grid we are working to retrieve the climatic variables for
grid_size <- '50'
grid_units <- 'km'

# Year period over which to work
start_year <- 1985
end_year <- 1993

# Filepath to datafiles with basename
# era_fp <- 'D:/carlota/projects/data/ERA5Land_cerrado/{year}{month_str}.{variable}_cerrado.tif'
# era_fp <- '../../data/raw/ERA5Land_cerrado/{year}{month_str}.{variable}_cerrado.tif'
era_fp <- '../../data/raw/ERA5Land_cerrado/{year}{month_str}.{variable}.tif'


# Output filepath
out_fp <- '../../data/processed/climate_data/era5land_grid{grid_size}{grid_units}_{year}.csv'


# ---------------- Reading grid ------------------------------------------------

# Filepath
grid_fp <- glue('../../data/shapes/cerrado_grid_{grid_size}{grid_units}_cropped.shp')
                
# Read the grid into memory
grid <- st_read(grid_fp)

# TODO: check grid columns here and discard useless

# ----------------- Setup cluster ----------------------------------------------
cl <- makeCluster(38)
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
    # grid_month <- precipitation(grid_month, era_fp, year, month)
    
    # 2. Temperature
    grid_month <- temperature(grid_month, era_fp, year, month)
    
    # 3. Relative Humidity
    grid_month <- RH(grid_month, era_fp, year, month)
    
    # 4. Vapour Pressure Deficit
    grid_month <- VPD(grid_month, era_fp, year, month)
    
    # 5. Potential Evapotranspiration
    grid_month <- potentialET(grid_month, era_fp, year, month)
    
    # 6. Actual Evapotranspiration
    grid_month <- actualET(grid_month, era_fp, year, month)
    
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