# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 11:02:25 2023

@author: scat8298
"""

import gdown
import pandas as pd

list_files = pd.read_csv('./MBlulc_download_c60_1986-2020.csv')


for file in list_files['files']:
    print('https://drive.google.com/uc?id=' + str(file))
    gdown.download(url= 'https://drive.google.com/uc?id=' + str(file))
