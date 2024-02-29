# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 18:52:52 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script containing two functions to calculate the overlap between two points 
of clouds by means of the Bhattacharyya distance.

"""

"""
Imports
"""
import numpy as np
import numpy.linalg as la
from sklearn.metrics.pairwise import euclidean_distances

"""
Functions
"""

def bhattacharyya_gaussian_distance(distribution1: "dict", distribution2: "dict",) -> int:
    """ Estimate Bhattacharyya Distance (between Gaussian Distributions)
    
    Args:
        distribution1: a sample gaussian distribution 1
        distribution2: a sample gaussian distribution 2
    
    Returns:
        Bhattacharyya distance
    """
    mean1 = distribution1["mean"]
    cov1 = distribution1["covariance"]

    mean2 = distribution2["mean"]
    cov2 = distribution2["covariance"]

    cov = (1 / 2) * (cov1 + cov2)
    
    
    
    # print('The problem is here?')
    # print(la.det(cov) / np.sqrt(la.det(cov1) * la.det(cov2)))
    
    # print(cov1)
    # print(la.det(cov1))
    # print(cov2)
    # print(la.det(cov2))
    # print(cov)
    # print(la.det(cov))

    T1 = (1 / 8) * (
        np.sqrt((mean1 - mean2) @ la.inv(cov) @ (mean1 - mean2).T)
    )
    T2 = (1 / 2) * np.log(
        la.det(cov) / np.sqrt(la.det(cov1) * la.det(cov2))
    )

    return T1 + T2


def overlap_BD(data1, data2):
    # This function calculates the Bhattarcchyya distance betweeen two sets of points
    
    # Calculate the empirical mean and covariances
    data_1 = data1.to_numpy()
    mean_1 = np.mean(data_1, axis = 0)
    cov_1 = np.cov(data_1, rowvar = False)

    
    data_2 = data2.to_numpy()
    mean_2 = np.mean(data_2, axis = 0)
    cov_2 = np.cov(data_2, rowvar = False)
    
    # Getting the diagonal elements (the sigmas)
    # sigmas_1 = np.diagonal(cov_1)
    # sigmas_2 = np.diagonal(cov_2)
    
    # Are they separable??
    # sep = (mean_1 - mean_2)**2/(sigmas_1**2 + sigmas_2**2)

    
    # Distribution 1
    distribution_1 = {
        "mean": mean_1,
        "covariance": cov_1,
    }
    
    
    # Distribution 2
    distribution_2 = {
        "mean": mean_2,
        "covariance": cov_2,
    }
    
    distance = bhattacharyya_gaussian_distance(distribution_1, distribution_2)
    
    # print('distnace' , distance)
    
    return distance

def overlap_1d(data1, data2):
    # Function that determines if two clusters should be merged or split based on
    # the overlap on the first PCA axes
    
    # Subset the relevant information
    data1 = data1.loc[:, 'PC1'].to_numpy().reshape(-1, 1)
    data2 = data2.loc[:, 'PC1'].to_numpy().reshape(-1, 1)
    
    # The merged cluster
    dataM = np.concatenate((data1, data2))
    
    # Calculate the maximum distance between points inside each cluster
    distM = np.max(euclidean_distances(dataM, dataM))
    dist1 = np.max(euclidean_distances(data1, data1))
    dist2 = np.max(euclidean_distances(data2, data2))
    
    merge = distM < dist1 + dist2
    split = not merge
    
    return merge, split
    
    
    