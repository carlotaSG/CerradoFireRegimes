# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 17:01:52 2023

@author: scat8298
"""

from multiprocessing import Pool

data = [1,2,3,4,5]

def dummy(d):
    return d+1, d+2
    

if __name__ == '__main__':
    
    with Pool(2) as p:
        output, outt = p.map(dummy, data)
    
    print(output)