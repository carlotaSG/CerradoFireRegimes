# var_exploration.R
#
# Created: 07/06/2023
# Author: Carlota Segura-Garcia
# Email: carlota.seguragarcia@ouce.ox.ac.uk
#
# Function that creates a bivariate plot with histograms of the variables in the
# diagonal (using the function panelutils() from "Numerical Ecology with R").
#
# ARGUMENTS
#
#

source('../../scripts/utils/panelutils.R')
source('../../scripts/utils/cleanplot.pca.R')

read_format_input <- function(data_fp, period){
  
  data <- read.csv(data_fp)
  
  rownames(data) <- data$id
  data$id <- NULL
  
  data$period <- NULL
  
  return(data)
}

calculate_log <- function(data, vars){
  
  for (var in vars){
    data[paste0( 'log_', var)] <- log(data[var])
  }
  
  return(data)
}

bivariatePlots_withLogs_fire <- function(data, log_vars = NULL, exc_vars = NULL, title_id = '') {
  # Function that creates a bivariate plot with histograms of the variables in the
  # diagonal (using the function panelutils() from "Numerical Ecology with R").
  
  # If there are any variables for which the user wants to calculate the logarithm
  if (length(log_vars) > 0){
    data <- calculate_log(data, log_vars)
  }
  
  # If there are any variables for which the user wants to calculate the logarithm
  if (length(exc_vars) > 0){
    data <- data[ , !names(data) %in% exc_vars]
  }

  # Generating the bivariate plots with all variables
  pairs(data,
        panel = panel.smooth,
        diag.panel = panel.hist,
        main = paste('Bivariate plots of fire variables for', title_id))
}


pca_screeplot_biplot <- function(data, log_vars = NULL, exc_vars = NULL, sites.labels = TRUE){
  
  # If there are any variables for which the user wants to calculate the logarithm
  if (length(log_vars) > 0){
    data <- calculate_log(data, log_vars)
  }
  
  # Calculate the PCA based on a correlation matrix (?)
  # Using scale = TRUE to standardise the variables
  data.pca <- rda(data[ , !names(data) %in% exc_vars], scale = TRUE)
  
  
  # Creating a scree plot to visualise the importance of each eigenvalue and  
  # in comparison to a broken stick model
  screeplot(data.pca, bstick = TRUE, npcs = length(data.pca$CA$eig))
  
  
  # Biplots of sites and variables
  par(mfrow = c(1,2))
  cleanplot.pca(data.pca, scaling = 1, mar.percent = 0.08, label.sites = sites.labels)
  cleanplot.pca(data.pca, scaling = 2, mar.percent = 0.5, label.sites = sites.labels)
}
