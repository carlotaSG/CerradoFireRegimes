# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 19:11:19 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that calculates the Bayesian Information Criteria for a set of points 
classified into different clusters.
"""

"""
Imports
"""
import numpy as np
from scipy.spatial import distance

"""
Functions
"""

def compute_bic(X, labels, centers, nclusters):
    """
    Computes the BIC metric for a given clusters

    Parameters:
    -----------------------------------------
    kmeans:  List of clustering object from scikit learn

    X     :  multidimension np array of data points

    Returns:
    -----------------------------------------
    BIC value
    
    
    Following: https://stats.stackexchange.com/questions/90769/using-bic-to-estimate-the-number-of-k-in-kmeans
    """
    
    # Number of clusters
    m = nclusters #kmeans.n_clusters
    
    
    

    #size of data set
    N, d = X.shape
    # print(N, d)
    
    const_term = 0.5 * m * np.log(N) * (d+1)

    #compute variance for all clusters beforehand
    if nclusters > 1:
        # size of the clusters
        n = np.bincount(labels)
        
        cl_var = (1.0 / (N - m) / d) * sum([sum(distance.cdist(X[np.where(labels == i)], [centers[i]], 'euclidean')**2) for i in range(m)])
        # print(cl_var)
        
        BIC = np.sum([n[i] * np.log(n[i]) -
               n[i] * np.log(N) -
             ((n[i] * d) / 2) * np.log(2*np.pi*cl_var) -
             ((n[i] - 1)*d / 2) for i in range(m)]) - const_term
        
    elif nclusters == 1:
        
        # Size of the clusters
        n = N
        
        cl_var = (1.0 / (N - m) / d) * sum(sum(distance.cdist(X, [centers], 'euclidean')**2))
        
        
        BIC = np.sum(n * np.log(n) -
               n * np.log(N) -
             ((n * d) / 2) * np.log(2*np.pi*cl_var) -
             ((n - 1)*d / 2)) - const_term

    

    return(BIC)   