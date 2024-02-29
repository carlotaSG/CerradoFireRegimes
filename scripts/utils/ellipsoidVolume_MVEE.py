# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 19:06:48 2023

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Scripts containing tha functions necessary to calculate the volume of the 
minimum volume ellipse of a cloud of points, as well as the comparison 
criteria betweeen two data points tand the union of them (following 
Lughofer and Sayed-Mouchaweh 2015).
"""

"""
Imports
"""
import numpy as np
import numpy.linalg as la

"""
Functions
"""

def mvee(points, tol = 0.001):
    """
    Find the minimum volume ellipse.
    Return A, c where the equation for the ellipse given in "center form" is
    (x-c).T * A * (x-c) = 1
    
    (04/08/2023)  https://stackoverflow.com/questions/14016898/port-matlab-bounding-ellipsoid-code-to-python
    """

    points = np.asmatrix(points)
    
    N, d = points.shape
    
    Q = np.column_stack((points, np.ones(N))).T
    
    err = tol+1.0
    u = np.ones(N)/N
    while err > tol:
        # assert u.sum() == 1 # invariant
        X = Q * np.diag(u) * Q.T
        M = np.diag(Q.T * la.inv(X) * Q)
        jdx = np.argmax(M)
        step_size = (M[jdx]-d-1.0)/((d+1)*(M[jdx]-1.0))
        new_u = (1-step_size)*u
        new_u[jdx] += step_size
        err = la.norm(new_u-u)
        u = new_u
    c = u*points
    A = la.inv(points.T*np.diag(u)*points - c.T*c)/d    
    
    # return np.asarray(A), np.squeeze(np.asarray(c))
    return A, c



def calculate_ellipsoidVolume(data, nk):
    # Data is alread a numpy array with the axis coordinates as columns
    # and the appropriate dimensionality if it had to be reduced
    
    # First, obtain the matrix of inequalities A, from A x = c
    A, _ = mvee(data)
    
    # print(type(A))
    
    # Calculate eigenvalues
    # eigenvalues, _ = la.eig(A.T*A) # THIS WAS AN ERROR!! FOR SOME REASON IT WAS YIELDING NEGATIVE EIGENVLAUES!!
    # eigenvalues, _ = la.eig(A)
    _, eigenvalues, _ = la.svd(A)
    
    # print(eigenvalues)
    
    # Factor rk
    rk = (1 - 1/(nk + 1))**(-4)
    
    # print(eigenvalues)
    # print('negative eigenvalues')
    # print(np.any(eigenvalues  < 0))
    # print(np.sqrt(rk/eigenvalues))
    # print(rk)
    
    # Volume of the ellipse (without factors)
    v = np.prod(np.sqrt(rk/eigenvalues))
    
    return v
    


def compare_ellispsoidVolumes(data_m, data_s1, data_s2, ndim, clust_ids = [1,2]):
    
    # Calculating the volume for the initial cluster
    
    # Convert to numpy array
    data_m = data_m.to_numpy()
    # The number of observations in cluster
    n = len(data_m)
    
    V_merge = calculate_ellipsoidVolume(data_m, n)
    
    
    # Calculating the volume for the first split cluster
    # Subsetting the relevant data
    data_s1 = data_s1.to_numpy()
    # The number of observations in cluster
    n = len(data_s1)
    
    V_prop1 = calculate_ellipsoidVolume(data_s1, n)
    
    
    # Calculating the volume for the first proposed cluster
    # Subsetting the relevant data
    data_s2 = data_s2.to_numpy()
    # The number of observations in cluster
    n = len(data_s2)
    
    V_prop2 = calculate_ellipsoidVolume(data_s2, n)
    
    
    inequality = ((ndim-2) * (V_prop1 + V_prop2)) / V_merge > 1.
    
    # Merge if inequality is true
    merge = inequality
    
    # Split if inequality is not true
    split = not merge
    
    return (merge, split)

def points_insideEllipsoid(points, cloud_points):
    # Function that obtains the matrix A (of inequalities) and the centre of an 
    # ellipsoid. Then, it returns for a list of points their "belonging" to the
    # ellipsoid
    #
    # Idea from (08/08/2023)
    # https://stackoverflow.com/questions/34376477/check-if-point-lies-inside-mulit-dimensional-ellipsoid
    #
    # A point is found within an ellipsoid if (x-c).trans A (x-c) < 1 + tolerance
    
    # tol_mvee = 0.01;
    # tol_dist = 0.1;

    A, c = mvee(cloud_points)
    

    # Initialising list to store the "belonging"
    d = []
    
    # For each point in list
    for i in range(len(points)):
        x   = points.iloc[i,:].transpose().to_numpy()
        x_c = x-c
        
        di = (x_c * (A * x_c.T))
        di = np.squeeze(np.array(di))
        d.append(di)
        # d.append(np.dot(x_c, np.dot(A, x_c.T))) 
        
    return d


def check_pointsInEllipsoid(list_data, list_labels):
    
    # First, determine which of the two has the smallest number of data points
    if len(list_data[0]) > len(list_data[1]):
        idmin = 1
        idmax = 0
    else:
        idmin = 0
        idmax = 1
        
    # Now we can calculate the index of "belonging" of each point in the small
    # cloud in the large cloud's MVEE

    out = points_insideEllipsoid(list_data[idmin], list_data[idmax])
    
    
    if sum([i <= 1.5 for i in out])/len(list_data[idmin]) >= 0.3:
        # In this case, we merge
        merge = True
        split = False
        
        # The index to assign is that of the larger cluster
        merge_label = list_labels[idmax]
        
    else:
        merge = False
        split = True
        merge_label = 0
        
    return merge, split, merge_label





