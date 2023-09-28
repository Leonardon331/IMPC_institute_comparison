# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 22:18:48 2023

@author: Leonard Borchardt
"""
### function not important


import os
import pandas as pd


#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')

#%%

def all_gene_controls():
    def all_gene_controls_list(path: str = r'./parameters', controls: list = []):
            
        contains = 'controls' # substring that the filename should contain
        
        for item in os.listdir(path):
            item_path = os.path.join(path, item)
            if os.path.isfile(item_path) and contains in item:
                df = pd.read_csv(item_path)
                controls.append(df)
                print(item, "appended")
                
            if os.path.isdir(item_path):
                all_gene_controls_list(item_path, controls)
        
        return controls
                
    controls = all_gene_controls_list()
    

    if len(controls) != 0:      
        controls = pd.concat(controls,ignore_index=True)
        controls = controls.sort_values(by=["parameter_stable_id", "sex","phenotyping_center"]).reset_index(drop=True)        
        controls.to_csv(r'./genes/all_gene_controls.csv', index=False)
    
#%%
all_gene_controls()
