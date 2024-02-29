# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 19:19:50 2021

@author: Carlota Segura-Garcia
@email: carlota.seguragarcia@ouce.ox.ac.uk

Script that retrieves the list of files in a directory. 

Optionally, it generates 
an "output" list of files with the same name as original files but with a user-provided
extension, and an output directory.
"""

"""
Imports
"""
import os

import csv


"""
Functions
"""

def fileList(inDir, outDir=None, suffix=None, in_extensions=None):
    """
    Function that takes a path to a directory as input and retrieves the list of files
    inside the given directory.
    
    Optionally, it returns a tuple with the above list of files as first element, 
    and a list of filapths with the output directory as path and the same filenames
    as in the input directory as second element.
    
    Optionally, if a suffix is passed, this suffix is added to the output filenames right before
    the file type.
    
    Optionally, if input extensions are passed, this function will only incpororate into the list 
    those filenames with the extension(s) in in_extensions.

    Parameters
    ----------
    inDir : string,
        The path to the directory where the files that we want to list are.
    outDir : string, optional
        The path to the folder where we want the output filenames to be written to. The default is None.
    suffix : string, optional
        A suffx we want to add to the output filenames. The default is None.
    in_extension: string, optional
        List of filename extensions that we want to read.

    Returns
    -------
    fp_list : list of string, or (list of strings, list of strings)
        List with the filepath of the files in the inDir directory. Optionally, if outDir is passed,
        it returns a tuple with the list of filepaths in the input directory, and the list of filepaths to 
        the output directory (outDir).

    """
    
    # If inDir has no '/' as last character, then add it
    if inDir[-1] != '/':
        inDir = inDir + '/'
    
    # Creating list of filepaths in directory inDir
    #   Encoding inDir
    encoded_inDir = os.fsencode(inDir)
    #   Getting encoded list of filenames
    encoded_infn = os.listdir(encoded_inDir)
    #   Decoding list of filenames to human readable
    decoded_infn = [os.fsdecode(filename) for filename in encoded_infn]
    
    #   If a list of file extensions is passed, drop those filenames with extensions not listed
    if in_extensions != None:
        decoded_infn = [file for file in decoded_infn if file.endswith(in_extensions)]

    #   Creating list of filepaths in directory inDir
    decoded_infp = [inDir + filename for filename in decoded_infn]
    

    
    # If an output directory is passed, then a list of output files is generated
    if outDir != None:
        # If outDir has no '/' as last character, then add it
        if outDir[-1] != '/':
            outDir = outDir + '/'
            
        # List of output files - same name as files in input directory inDir, but outDir as filepath
        decoded_outfp = [outDir + filename for filename in decoded_infn]
        
        # If suffix provided, adding suffix right before the file extension
        if suffix != None:
            decoded_outfp = [filepath[:-4] + suffix + filepath[-4:] for filepath in decoded_outfp]
    
    
    # List of filepaths to return
    if outDir == None:
        # If no out directory, only the files in the input are returned
        fp_list = decoded_infp
    elif outDir != None:
        # If out directory, then return tuple with list of filepaths in input directory and list of output filepaths
        fp_list = (decoded_infp, decoded_outfp)
        
    return fp_list


def createListOutputFiles(in_fp, outDir, suffix=None, out_extension=None):
    # If outDir has no '/' as last character, then add it
    if outDir[-1] != '/':
        outDir = outDir + '/'
        
    # If output directory does not exist, print warning and create it
    if not os.path.exists(outDir):
        print('Warning: output directory did not exist, creating it.')
        print('Output directory: {}'.format(outDir))
        
        os.makedirs(outDir)
        
    # List of output files - same name as files in input directory inDir, but outDir as filepath
    out_fp = [outDir + os.path.basename(filepath) for filepath in in_fp]
    
    # If suffix provided, adding suffix right before the file extension
    if suffix != None:
        out_fp = [filepath[:-4] + suffix + filepath[-4:] for filepath in out_fp]
        
    
    # If extensions provided, substituting extension
    if out_extension != None:
        out_fp = [filepath[:-4] + out_extension for filepath in out_fp]
        
    return out_fp



def readListFiles_fromCSV(csv_fp):
    """
    Function that gets a csv - presumably containing fileapths or list of years - 
    and returns a list of filepaths if csv contains one column, or a list of
    tuples - one for each row - if csv contians two columns.

    Parameters
    ----------
    csv_fp : string
        Filepath to the CSV file containing the information.

    Returns
    -------
    list_fp : list of strings or list of tuples
        List of strings if csv contains one columns, list of tuples if csv contains
        two columns. One tuple per row.

    """
    
    #fileList_csv = './server_orders/MBfogo_coordinates_node9.csv'
    # UPDATE (16/03/2023)
    with open(csv_fp, newline='') as f:
        reader = csv.reader(f)
        
        # Counting number of columns
        ncol = len(next(reader))
        f.seek(0)
        
        # If reader only has one column, return a list
        if ncol == 1:
            list_fp = [element[0] for element in reader]
        if ncol == 2:
            list_fp = [(element[0], element[1]) for element in reader]   
    
    # Old code - for legacy purposes
    # with open(csv_fp, newline='') as f:
    #     reader = csv.reader(f)
    #     list_fp = [element[0] for element in reader]
        
    return list_fp

def readListFiles_withConditions(inDir, csv_fn_conditions):
    
    # Reading list of string conditions from csv file
    with open(csv_fn_conditions, newline='') as f:
        reader = csv.reader(f)
        conditions = [str(element[0]) for element in reader]
        
    # Reading all files in directory
    list_fp = fileList(inDir)
    
    # Keep only those files that contain one of the conditions:
    list_fp = [file_fp for condition in conditions for file_fp in list_fp if condition in file_fp]
    
    # Return list of filepaths that meet one of the conditions
    return list_fp
            
def selecting_fp_withCondition(list_fp, condition_list):
    """
    Function that returns all the filepaths in a list that contain at least one 
    string in a list of string conditions.

    Parameters
    ----------
    list_fp : list of strings
        List of filepaths. Those filepaths containing at least one of the strings
        in condition_list will be kept.
    condition_list : list of strings
        List of conditions (in string) format used to subset list_fp.

    Returns
    -------
    list_fp : list of strings
        List of filepaths where each filepath contains at least one of the strings
        in condition_list.

    """
    # Subsetting list
    list_fp = [fp for condition in condition_list for fp in list_fp if condition in fp]
    
    # Keeping only unique entries
    list_fp = list(set(list_fp))
    
    return list_fp      

def listAllFiles_dirAndSubdir(filepath, include_ext=None):
    # Function that gets the fileapths to all files in directory including subdirectories. 
    # Optionally it excludes files with extensions in exclude_ext
    
    list_fp = []
    for path, subdirs, files in os.walk(filepath):
        for name in files:
            # If file extension requested and name contains extension, then append 
            # to list
            if not (include_ext == None) and (include_ext in name):
                list_fp.append(os.path.join(path, name))
    
    return list_fp
    
            
