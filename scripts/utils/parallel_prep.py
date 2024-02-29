# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 14:37:43 2022

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that helps preparing the input to a parallel computation. 

The following functions can be called form other scripts:
    - list_slices() - Function that, given a total number of rows and a number of 
                      partitions, returns a list of slices of indices corresponding 
                      to sets of rows in a dataframe. This allows calling script 
                      to operate in parallel over these sets of rows.
"""

"""
Imports
"""
import numpy as np

"""
Functions
"""
def list_slices(ntotal, partitions):
    """
    Function that, given the total number of rows in a file and a number of partitions,
    returns the slices of the range (0, ntotal), as many slices as the number of partitions.
    
    All slices have the same length except the last one, which includes the
    mode (the residual) in its length.   

    Parameters
    ----------
    ntotal : integer
        The numer of rows in the data (the range of values).
    partitions : integer
        The number of partitions in which to divide the rows (range of values).

    Returns
    -------
    slices : list of tuples
        List of slices of indices (row names) to be used when partitioning the 
        data for parallel computation.

    """
    
    # Number of elements per slice rounded down
    nelements = int(ntotal/partitions)
    
    # The mode or residual
    the_mod = ntotal % partitions
    
    # The right hand side of the slice
    right = np.arange(0, nelements * partitions, nelements)
    
    # The left hand side is the same as the right but starting from the second element
    left = right[1:len(right)]
    # And adding the very last element
    left = np.append(left, left[len(left)-1] + nelements + the_mod)
    
    # Finally, the slices as a list of tuples
    slices = list(zip(right, left))
    
    return slices

