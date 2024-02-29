# read_format_grid.R
#
# Created: 09/06/2023
# Author: Carlota Segura-Garcia
# Email: carlota.seguragarcia@ouce.ox.ac.uk
#
# Function that creates a bivariate plot with histograms of the variables in the
# diagonal (using the function panelutils() from "Numerical Ecology with R").
#
# ARGUMENTS
#
#
library(sf)

get_grid <- function(grid_id, area_thresh) {
  
  # Reading the grid data
  grid <- st_read(glue('../../data/shapes/cerrado_grid_{grid_id}_cropped.shp'),
                  quiet = TRUE)
  
  # Selecting relevant columns
  grid <- grid %>% 
    select(id, X.area_crop, geometry)
  
  
  # Selecting only those cells with more than x% of their area in the Cerrado
  grid <- grid[((grid$X.area_crop >= area_thresh) | (grid$id %in% c(567, 568))) & !(grid$id %in% c(1451)), ]
  
  return(grid)
}