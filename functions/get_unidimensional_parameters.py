# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 10:28:30 2023

@author: Leonard Borchardt
"""

import os
import warnings
import pandas as pd
import numpy as np

#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')
#%%

def get_unidimensional_parameters(gene_symbols: list):
    '''
    checks given parameters for unidimensionality and returns a list with unidimensional parameters only,
    i.e., no series or categorical (and so on) parameters

    import modules:
        import os
        import warnings
        import pandas as pd

    Parameters
    ----------
    gene_symbols : list of str, gene symbols

    Returns
    -------
    None
    
    automatically assigns unidimensional-parameter-only list to global variable
    saves the list as a csv data frame, too

    '''
    
    # ignore warnings
    warnings.filterwarnings("ignore")
    
    for gene_symbol in gene_symbols:
        print(gene_symbol)
        read = f"./genes/{gene_symbol}/{gene_symbol}_most_center_ps_ids.csv"
        save = f"./genes/{gene_symbol}/{gene_symbol}_most_center_ps_ids_unidimensional.csv"
        
        parameter_stable_ids = pd.read_csv(read).parameter_stable_id.tolist()
        
        df = pd.read_csv(r'./genes_all_columns/all_genes_all_columns.csv')
        
        #choosing parameter data for parameters that are part of the parameter_stable_ids list and unidimensional
        df = df[(df['parameter_stable_id'].isin(parameter_stable_ids))&(df.observation_type == "unidimensional")].parameter_stable_id.unique()
        
        #remove IMPC_BWT_001_001 (body weight) from array (not so easy to compare between centers)
        df = np.delete(df, np.where(df == "IMPC_BWT_001_001"))
        
        #array to frame
        df = pd.DataFrame(df, columns = ["parameter_stable_id"])
        
        df = df.sort_values(by=["parameter_stable_id"]).reset_index(drop = True)
        
        df.to_csv(save,index=False)
        globals()[f"{gene_symbol}_most_center_ps_ids_unidimensional"] = df
    
#%%

get_unidimensional_parameters(["Prkab1", "Dbn1", "Ap4e1", "Nxn"])