# cluster_exploration.R
#
# Created: 07/06/2023
# Author: Carlota Segura-Garcia
# Email: carlota.seguragarcia@ouce.ox.ac.uk
#
#
# ARGUMENTS
#
#

library(gclus)
library(tidyr)
library(cowplot)


plot_fusion_levels <- function(clust_obj, max_k, title = 'Fusion levels'){

  # Graph of fusion level values
  plot(tail(clust_obj$height, max_k - 1),
       max_k:2,
       type = 'S',
       main = title,
       ylab = 'k (number of clusters)',
       xlab = 'h (node height)',
       col = 'grey')
  
  # Adding the labels indicating the number of clusters
  text(tail(clust_obj$height, max_k - 1),
       max_k:2,
       max_k:2,
       col = 'red',
       cex = 0.8)
}


plot_optimal_silhouette <- function(clust_obj, dis_matrix, max_k,
                                    title = 'Silhouette-optimal number of clusters'){
  
  # Creating 0 vector for each number of clusters solution to store the 
  # silhouette value for each solution
  sil_k <- numeric(max_k)
  
  for (k in 2:max_k){
    # Calculate the silhoutte width of each site for this k solution
    sil <- silhouette(cutree(clust_obj, k = k), dis_matrix)
    
    # Summarising to have the average silhouette value per cluster, and then getting
    # the average silhouette among the clusters
    sil_k[k] <- summary(sil)$avg.width
  }
  
  # We consider the best solution to be the one with largest average silhouette width
  k_best <- which.max(sil_k)
  
  # Plotting the average silhouette width per solution
  plot(1:max_k,
       sil_k,
       type = 'h',
       main = title,
       xlab = 'k (number of clusters)',
       ylab = 'Average silhouette width')
  
  # Highlighting and labeling the optimal number of clusters according to the 
  # average silhouette width
  axis(1,
       k_best,
       paste('optimum', k_best, sep = '\n'),
       col = 'red',
       font = 2,
       col.axis = 'red')
  points(k_best,
         max(sil_k),
         pch = 16, col = 'red', cex = 1.5)
}


map_k_clusters <- function(clust_k_classification, k, clust_id_col, clust_col, grid, pal, plt_contour = FALSE){
  
  
  clust_k <- clust_k_classification[ , c(clust_id_col, clust_col)]
  
  colnames(clust_k) <- c('id', 'cluster')
  
  
  # Adding data to grid 
  grid <- grid %>% 
    left_join(clust_k,
              by = 'id')
  
  grid <- grid %>% 
    drop_na()
  
  
  # Creating the map
  #-------------------
  
  # Creating the labels for map's legend
  labs_plot <- levels(factor(grid$cluster))#1:k
  
  # Palette
  # pal <- hcl.colors(k+1, "Inferno", rev = TRUE, alpha = 0.7)
  
  
  # Now we could map the clustering solution
  map <- ggplot() + 
    # Add choropleth overlay
    geom_sf(data = grid,
            aes(fill = factor(cluster)), color = NA) +
    # Labs
    labs(#title = glue("k = {k}"),
         # subtitle = "Hierarchical clustering, Euclidean distance, normalised",
         # caption = "Cluster ID",
         fill = "") +
    # Custom palette
    scale_fill_manual(values = pal,
                      drop = FALSE,
                      na.value = "grey80",
                      label = labs_plot,
                      # Legend
                      guide = guide_legend(direction = "horizontal",
                                           nrow = 1,
                                           label.position = "bottom")) +
    # Theme
    theme_void() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 9),
          legend.key.width = unit(3, "mm"),
          legend.key.height = unit(3, "mm"))
  
  if (plt_contour == TRUE){
    map <- map +
      # theme(panel.background = element_rect(fill = NULL, colour = 'red', linewidth = 2))
      theme(plot.title = element_text(color = 'red'))
  }
  
  
  return(map)
}


map_various_clusters <- function(clust_obj, clusters, k_opt, grid,
                                 n_col = 2, n_row = 2,
                     title = 'Maps of different clusters',
                     subtitle = ''){
  
  plots <- list()
  
  for (i in 1:length(clusters)){
    
    if (clusters[i] == k_opt){
      plot_contour = TRUE
    } else {
      plot_contour = FALSE
    }
    
    # Create map for k clusters solution
    map_k <- map_k_clusters(clust_obj, clusters[i], grid, plot_contour)
    
    plots[[i]] <- map_k
    
    
  }
  

  themaps <- plot_grid(plotlist = plots, ncol = n_col, nrow = n_row)
  
  thetitle <- ggdraw() + 
    draw_label(title, size = 18, x = 0, hjust = 0) + 
    draw_label(subtitle, size = 8,
               y = 0.15, x = 0, hjust = 0) 
  
  plot_grid(thetitle, themaps, ncol = 1, rel_heights = c(0.15, 0.85))+ 
    theme(plot.margin = margin(0, 0, 0, 7),
          panel.background = element_rect(fill = 'white', colour = 'white'))
}