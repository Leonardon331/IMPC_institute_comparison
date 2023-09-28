# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 08:36:05 2023

@author: Leonard Borchardt
"""

import os
import warnings
import pandas as pd

#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')

#%%

def concat_csv_files(files_name: str, files_path: str = r'./', save_path: str = r'./', save_name: str = 'concatenated_csv'):
    '''
    concatenates specific csv files
    
    import modules:
        import os
        import warnings
        import pandas as pd

    Parameters
    ----------
    files_name: substring that the file names should contain
    files_path: path that incl. all files, default: cwd, optional
    save_path: path for saving the concatenated data frame, default: cwd, optional
    save_name: name of the concatenated data frame for saving, default: 'concatenated_csv', optional
    
    Returns
    -------
    None
    
    '''
    
    # ignore warnings
    warnings.filterwarnings("ignore")
    
    if not files_path.endswith('/'):
        files_path += r'/'
        
    if not save_path.endswith('/'):
        save_path += r'/'
    
    if not save_name.endswith('.csv'):
        save_name += r'.csv'
    
    def all_files_list(files_name: str, files_path: str, files: list = []):
        for item in os.listdir(files_path):
            item_path = os.path.join(files_path, item)
            if os.path.isfile(item_path) and files_name in item:
                df = pd.read_csv(item_path)
                files.append(df)
                print(item, "appended")
                
            if os.path.isdir(item_path):
                all_files_list(files_name, item_path, files)
        
        return files
                
    files = all_files_list(files_name, files_path)
    

    if len(files) != 0:      
        files = pd.concat(files,ignore_index=True)
        files = files.sort_values(by=["gene_symbol","parameter_stable_id","zygosity","sex","phenotyping_center"]).reset_index(drop=True)
        files.to_csv(save_path+save_name, index=False)
        print("all data frames concatenated")
    
    else:
        print("no concatenation")
        print("check files_name")
    
#%%
concat_csv_files(files_name = '_category_proportions', save_path = r'./overall_analysis/categorical_data', save_name = 'category_proportions')