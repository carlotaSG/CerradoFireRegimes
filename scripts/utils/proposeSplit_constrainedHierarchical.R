#...............................................................................
#                                                                              .
#  Created: 04/08/2023                                                         .
#  Author: Carlota Segura-Garcia                                               .
#  Email: carlota.seguragarcia@ouce.ox.ac.uk                                   .
#                                                                              .
#  Script that applies a Spatially Constrained Hierarchical Algorithm          .
#  (following Guenard and Legendre 2022) to create 2 clusters. Returns a       .
#  dataframe with the classification results for the two solutions. In the     .
#  output dataframe, cells are identified by their grid id and by the ordered  .
#  (1:length(grid)) id used by the clustering functionthe assigned cluster,    . 
#  either 1 or 2.                                                              . 
#                                                                              .
#  This script is to be run by calling it from a Python script and it requires .
#  two arguments: the original cluster label fo the data points to be split,   .
#  and the period in which we are working.                                     .
#                                                                              .
#...............................................................................

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                  Imports                                 ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(glue)
library(ade4)
library(adespatial)
library(dplyr)
library(magrittr)
library(vegan)
library(cluster)
library(tidyr)
library(clValid)



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                  Functions                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


read_format_input <- function(data_fp, period){
  
  data <- read.csv(data_fp)
  
  rownames(data) <- data$id
  data$id <- NULL
  
  data$period <- NULL
  
  return(data)
}

id_translation <- function(original_id){
  
  original_to_ordered <- setNames(1:length(original_id), original_id)
  ordered_to_original <- setNames(original_id, 1:length(original_id))
  
  return(list(original_to_ordered = original_to_ordered, ordered_to_original = ordered_to_original))
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                    MAIN                                  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#..........................User inputs...........................

args <- commandArgs(trailingOnly = TRUE)

# Periods to work over
period <- args[1]

# The cluster that we are working on
cluster <- args[2]

# The number of years per period
nyears <- args[3]

# # Whether to standardise the PCA scores again or not
# timesStandard <- args[3]

# Filepath to edge list
edge_list_fp <- glue("../../analysis/{nyears}_year_periods/temp/period{period}_cluster{cluster}_edges-list.csv")

# Filepath to PCA scores data
pca_scores_fp <- glue("../../analysis/{nyears}_year_periods/temp/period{period}_cluster{cluster}_toSplit.csv")

# Output file to store the data
output_fp <- glue("../../analysis/{nyears}_year_periods/temp/period{period}_cluster{cluster}_split.csv")

#.......................Clustering process.......................

# First, read the edge list (it is the same for all periods)
edge_list <- read.csv(edge_list_fp)

# Read the data 
data <- read_format_input(pca_scores_fp, p)


# To use the clustering algorithm, we have to rename the cells id to start fomr 1
translation <-  id_translation(rownames(data))

# Storing the two ways of translation
original_to_ordered <- translation[[1]]
ordered_to_original <- translation[[2]]

edge_list$from <- recode(edge_list$from, !!!original_to_ordered)
edge_list$to <- recode(edge_list$to, !!!original_to_ordered)

# Transform the id's in the data
rownames(data) <- recode(as.integer(rownames(data)), !!!original_to_ordered)


#................Fitting the clustering algorithm................



# Whether we standardise the PCA scores or not
timesStandard <- 'single'
if (timesStandard == 'double'){
  # Variable standardisation (using vegan package)
  data_st <- decostand(data, 'standardize', MARGIN = 2)
  
}else if (timesStandard == 'single'){
  data_st <- data
}


# Dissimilarity matrix using Euclidean distance (using adespatial package)
dis_mtx <- vegdist(data_st, 'euclidean')


# Spatially constrained clustering algorithm using the edge list
# (output from the clustering algorithm)
clust_out <- constr.hclust(d = dis_mtx,
                           links = edge_list[, c(1,2)])


# For k=2, obtain the classification
clust_out <- cutree(clust_out, k = 2)

# Transforming into a data frame for ease of use
clust_out <- data.frame(clust_id = as.integer(names(clust_out)), cluster = clust_out, row.names = NULL)

# Changing the cluster name
names(clust_out)[names(clust_out) == 'cluster'] <- glue('prop_clust')


# Add column with the original id (translating cell labels back)
clust_out$id <- recode(clust_out$clust_id, !!!ordered_to_original)

# Moving column to first position
clust_out %<>% select('id', everything())

# Dropping unnecessary column
clust_out <- clust_out %>% 
  select(-c('clust_id'))


# Saving to file
write.csv(clust_out, output_fp, row.names = FALSE)

