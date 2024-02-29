#...............................................................................
#                                                                              .
#  Created: 21/07/2023                                                         .
#  Updated: 26/10/2023                                                         .
#  Author: Carlota Segura-Garcia                                               .
#  Email: carlota.seguragarcia@ouce.ox.ac.uk                                   .
#                                                                              .
#  Script that contains the model specifications to calculate the PCA scores   .
#  from the transformed input variables, and subsequently run the clustering   .
#  algorithm.                                                                  .
#                                                                              .
#...............................................................................


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                            Model specifications                          ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Grid to work with
grid_id <- '30km'

# The variables to include in the analsyis (the non-transformed versions)
input_vars <- c('size_q99', 'size_q50', 
                'burned_area_n', 'number_fires_n', 'frequency_q99s')#'season_q50s', 


# The variables to include in the analysis (the transformed versions)
input_vars_transformed <- c('log_size_q99', 'log_size_q50', 
                            'log_burned_area_n', 'log_number_fires_n', 'log1_frequency_q99s')#'season_q50s', 

# Number of variables
n_vars <- length(input_vars)

# Number of years per period
nyears <- 9


# # Whether to work with missing data imputation and what number of fires       UPDATE 24/08/2023
# # to consider as "missing data"
# imputation_label <- 'No'                                                        # Options: No, 30, 40, 50


# Number of axes to cluster with (only relevant when running the clustering 
# algorithm)
n_PCAaxes <- 3 # n_vars


# # Whether to standardise the PCA scores before running the clustering algorithm  UPDATE 24/08/2023
# # (the transformed variable are always standardised BEFORE running the PCA. Here
# # we are deciding whether to standardise the PCA scores when calculating the 
# # dissimilarity matrix)
# timesStandard <- 'single'                                                       # Options: single, double


# Whether we are using a spatially constrained or an unconstrained algorithm
clustering_method <- 'const'                                              # Options: conStrained, unconstrained

# Connectivity
# Which edge list we consider using (whether 4 or 8 neighbors)
connectivity <- 8                                                               # Options: 8, 4



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                          Input/output data files                         ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ Data and PCA score files  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Directory with input data, PCA scores, and edge-list file
input_data_dir <- '../../data_new_grids/processed/clustering_algorithm/'


# Filepath to fire characteristics data
input_data_fp <- paste0(input_data_dir, '{nyears}_year_periods/input_vars_grid{grid_id}_period{p}.csv')

# Filepath to transformed fire characteristics data
transformed_input_data_fp <- paste0(input_data_dir, '{nyears}_year_periods/input_vars_transformed_grid{grid_id}_period{p}.csv')

# Filepath to calculated PCA scores
pca_scores_fp <- paste0(input_data_dir, '{nyears}_year_periods/fS_int_noSeas/PCAscores_scaling1_sites_input_vars_transformed_grid{grid_id}_period{p}.csv')

# Filepath to calculated PCA principal axes in feature space
pca_scores_species_fp <- paste0(input_data_dir, '{nyears}_year_periods/fS_int_noSeas/PCAscores_scaling1_species_input_vars_transformed_grid{grid_id}_period{p}.csv')

# File summarising the explanatory power of the PCA axes
pca_power_fp <- paste0(input_data_dir, 'PCAexpPower_scaling1_sites_input_vars_transformed_grid{grid_id}.csv')




##~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ PCA output plots  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~


# Plots that are exclusively output from the PCA scores fitting
PCA_output_dir <- '../../analysis_new_grids/{nyears}_year_periods/fS_int_noSeas/plot_pca/'

dir.create(glue(PCA_output_dir), showWarnings = FALSE)

# Biplot
# of scaling 1 (interpeting distances between data points)
biplot_scaling1_fp <- paste0(PCA_output_dir, 'biplot_scaling1_input_vars_transformed_grid{grid_id}_period{p}.png')

# of scaling 2 (interpeting distances between data points)
biplot_scaling2_fp <- paste0(PCA_output_dir, 'biplot_scaling2_input_vars_transformed_grid{grid_id}_period{p}.png')

# Screeplot
screeplot_fp <- paste0(PCA_output_dir, 'screeplot_input_vars_transformed_grid{grid_id}_period{p}.png')




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ Clustering algorithm  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Filepath to list of edges
edge_list_fp <- glue(paste0(input_data_dir, 'edge_list_connectivity{connectivity}_grid_{grid_id}_epsg5880.csv'))

# Identifier for the clustering algorithm method and specifications
# clustering_specs <- 'PCAscores_sECDF_{clustering_method}_{n_vars}vars_{n_PCAaxes}axes_imp{imputation_label}_{timesStandard}Stand_conn{connectivity}'

# Directory to store the clustering output
clustering_output_dir <- glue('../../analysis_new_grids/{nyears}_year_periods/fS_int_noSeas/')
clustering_output_plots_dir <- paste0(clustering_output_dir, 'plots_numClusters_grid{grid_id}/')

dir.create(glue(clustering_output_plots_dir), showWarnings = FALSE)

# File containing the classification for the various solutions
clustering_classification_fp <- paste0(clustering_output_dir, '{clustering_method}ClusteringClass_ksolutions_grid{grid_id}_period{p}.csv')


# The various plots analysing the solutions
dendrogram_fp <- paste0(clustering_output_plots_dir, 'dendrogram_grid{grid_id}_period{p}.png')
sil_fp <- paste0(clustering_output_plots_dir, 'silhouette_grid{grid_id}_period{p}.png')
CHC_WSS_fp <- paste0(clustering_output_plots_dir, 'calinski-harabasz_grid{grid_id}_period{p}.png')
fusionLevel_fp <- paste0(clustering_output_plots_dir, '_fusionLevel_grid{grid_id}_period{p}.png')
dunnIndex_fp <- paste0(clustering_output_plots_dir, 'dunnIndex_grid{grid_id}_period{p}.png')
bic_fp <- paste0(clustering_output_plots_dir, 'BIC_grid{grid_id}_period{p}.png')


##~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ Map and boxplot  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~

clustering_map_dir <- paste0(clustering_output_dir, 'plots_kSolutions_grid{grid_id}/')
dir.create(glue(clustering_map_dir), showWarnings = FALSE)
map_boxplot_fp <- paste0(clustering_map_dir, 'map_boxplot_grid{grid_id}_period{p}_clusters_{k}.png')
