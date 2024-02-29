#...............................................................................
#                                                                              .
#  Created: 03/11/2023                                                         .
#  Author: Carlota Segura-Garcia                                               .
#  Email: carlota.seguragarcia@ouce.ox.ac.uk                                   .
#                                                                              .
#  Script that reads all the explanatory variables, and returns all the data   . 
#  in one dataframe. Optionally, it includes the clustering classification for . 
#  the same period for a selection of solutions.                               .
#                                                                              .
#...............................................................................

read_explanatoryVariables <- function(nyears, grid_id, p_label, list_vars, clust = FALSE){
  
  #..................Reading and formatting data...................
  
  # Variables to keep from dataframe
  vars_to_keep <- append(c('id', 'period'), list_vars)

  # Explanatory variables
  # Land-use area and change
  lulc_dat <- read.csv(glue('../../data_new_grids/processed/summary_tables/{nyears}_year_periods/lulc_percArea_grid{grid_id}_allPeriods.csv'))
  lulc_dat <- lulc_dat[ , colnames(lulc_dat) %in% vars_to_keep]
  
  # Population density
  pop_dat <- read.csv(glue('../../data_new_grids/processed/summary_tables/{nyears}_year_periods/population_density_grid{grid_id}_allPeriods.csv'))
  pop_dat <- pop_dat %>% 
    rename(popdens = areaw_popdensity)
  pop_dat <- pop_dat[ , colnames(pop_dat) %in% vars_to_keep]
  
  # Livestock density
  lstock <- read.csv(glue('../../data_new_grids/processed/summary_tables/{nyears}_year_periods/livestock_density_grid{grid_id}_allPeriods.csv'))
  lstock <- lstock %>% 
    rename(lstkdens = areaw_dens)
  lstock <- lstock[ , colnames(lstock) %in% vars_to_keep]
  
  # Fragmentation
  frags <- read.csv(glue('../../data_new_grids/processed/summary_tables/{nyears}_year_periods/landscape_metrics_grid{grid_id}_allPeriods.csv'))
  frags <- frags[ , colnames(frags) %in% vars_to_keep]
  
  # Fragmentation at cell level
  cellfrags <- read.csv(glue('../../data_new_grids/processed/summary_tables/{nyears}_year_periods/landscape_cell_metrics_grid{grid_id}_allPeriods.csv'))
  # frags <- frags[c('id', 'period', 'patch_density', 'area_mn')]
  # print(cellfrags)
  cellfrags <- cellfrags[cellfrags['class_val'] == 'natural', ]
  cellfrags['class_val'] <- NULL
  colnames(cellfrags)[-c(1,2)] <- paste("cell", colnames(cellfrags)[-c(1,2)], sep = "_")
  cellfrags <- cellfrags[ , colnames(cellfrags) %in% vars_to_keep]
  
  # Explanatory variables
  # Climatic data
  clim_dat <- read.csv(glue('../../data_new_grids/processed/summary_tables/{nyears}_year_periods/clim_data_period_allVariables_grid{grid_id}_{nperiods}periods.csv'))
  clim_dat <- clim_dat[ , colnames(clim_dat) %in% vars_to_keep]
  
  
  # Merging into a single data frame
  data <- inner_join(lulc_dat, pop_dat, by = c('id', 'period'))
  data <- inner_join(data, lstock, by = c('id', 'period'))
  data <- inner_join(data, frags, by = c('id', 'period'))
  data <- inner_join(data, cellfrags, by = c('id', 'period'))
  data <- inner_join(data, clim_dat, by = c('id', 'period'))
  data <- data[data['period'] == p_label, ]
  data['period'] <- NULL

  # Topographic data
  topo_dat <- read.csv(glue('../../data_new_grids/processed/summary_tables/topography_allVars_grid{grid_id}.csv'))
  topo_dat <- topo_dat[ , colnames(topo_dat) %in% append(c('id'), list_vars)]
  data <- inner_join(data, topo_dat, by = 'id')
  
  if (clust == TRUE){
    # Cluster classification data
    clust_class <- read.csv(glue('../../analysis_new_grids/{nyears}_year_periods/fireSize_intersecting/constClusteringClass_ksolutions_grid{grid_id}_period{p}.csv'))
    
    # Adding to dataset
    data <- inner_join(clust_class, data, by = 'id')
  }
  
  return(data)
  
}


transform_explanatoryVariables <- function(data, list_vars, list_vars_trans){
  
  for (i in 1:length(list_vars)){
    if (list_vars_trans[i] == 'log'){
      data[paste0('log_', list_vars[i])] <- log(data[list_vars[i]])
      
    } else if (list_vars_trans[i] == 'sqrt'){
      data[paste0('sqrt_', list_vars[i])] <- sqrt(data[list_vars[i]])
      
    } else if (list_vars_trans[i] == 'log_rev'){
      data[paste0('log_', list_vars[i])] <- log(data[list_vars[i]] + abs(min(data[list_vars[i]])) + 1)
      
    } else if (list_vars_trans[i] == 'log1'){
      data[paste0('log1_', list_vars[i])] <- log(1 + data[list_vars[i]])
    }
  }
  
  return(data)
}

