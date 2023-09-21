# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 23:01:40 2023

@author: Leonard Borchardt
"""

import os

#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')
#%%

def rename_files(contains: str, change: str, contains_not: str = "     ", path: str = r'./'):
    '''
    deletes files starting from a given main directory that contain a certain string in their filename
    
    import os
    
    Parameters
    ----------
    path: str, the main directory where the function should search for files, default: cwd
    contains: str, string that the filename contains
    change: str, string that should replace contains
    contains_not: str, string that the filename not contains, optional

    Returns
    -------
    None.

    '''
    
    for item in os.listdir(path):
        item_path = os.path.join(path, item)
        if os.path.isfile(item_path) and contains in item and contains_not not in item:
            new_path = item_path.replace(contains, change)
            os.rename(item_path, new_path)
            print(item, "renamed")
            
        if os.path.isdir(item_path):      
            rename_files(contains, change, contains_not, item_path)

#%%
path = r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis\procedures'
rename_files("center_statistics", "ANOVA", path)