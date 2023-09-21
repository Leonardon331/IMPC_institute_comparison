# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 23:44:25 2023

@author: Leonard Borchardt
"""

import os

#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')
#%%

def delete_files(contains: str, contains_not: str = "     ", path: str = r'./'):
    '''
    deletes files starting from a given main directory that contain a certain string in their filename

    import os
    
    Parameters
    ----------
    contains : str, string that the filename contains
    contains_not: str, string that the filename not contains, optional
    path: str, the main directory where the function should search for files, default: cwd, optional

    Returns
    -------
    None.

    '''
    
    for item in os.listdir(path):
        item_path = os.path.join(path, item)
        if os.path.isfile(item_path) and contains in item and contains_not not in item:
            os.remove(item_path)
            print(item, "removed")
            
        if os.path.isdir(item_path):      
            delete_files(contains,contains_not,item_path)

#%%
delete_files(contains="kwr", contains_not="kwr_", path=r'./parameters')