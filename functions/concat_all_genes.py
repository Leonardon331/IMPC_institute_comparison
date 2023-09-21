# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 15:18:11 2023

@author: Leonard Borchardt
"""

#concatenate all gene_all_columns dfs

import os
import warnings
import pandas as pd



#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')

#%%

def concat_all_genes():
    '''
    concatenates all gene dfs with all columns
    
        run download_genes.py first
        
        import os
        import warnings
        import pandas as pd

    '''
    
    # ignore warnings
    warnings.filterwarnings("ignore")
    
    def all_genes_list(path: str = r'./genes_all_columns', genes: list = []):
            
        contains = 'all_columns' # substring that the filename should contain
        
        for item in os.listdir(path):
            item_path = os.path.join(path, item)
            if os.path.isfile(item_path) and contains in item:
                if item == "all_genes_all_columns.csv":
                    continue
                df = pd.read_csv(item_path)
                genes.append(df)
                print(item, "appended")
                
            if os.path.isdir(item_path):
                all_genes_list(item_path, genes)
        
        return genes
                
    genes = all_genes_list()
    

    if len(genes) != 0:      
        genes = pd.concat(genes,ignore_index=True)
        genes["date_of_birth"]= pd.to_datetime(genes["date_of_birth"])
        genes = genes.sort_values(by=["date_of_birth"], ascending = False).reset_index(drop=True)
        genes = genes.sort_values(by=["gene_symbol",,"procedure_name","parameter_stable_id","sex","phenotyping_center"]).reset_index(drop=True)
        genes.to_csv(r'./genes_all_columns/all_genes_all_columns.csv', index=False)
        print("all data frames concatenated")
    
    else:
        print("download genes first")
    
#%%
concat_all_genes()