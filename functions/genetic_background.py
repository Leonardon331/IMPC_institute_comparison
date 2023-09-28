# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 08:28:56 2023

@author: Leonard Borchardt
"""
#check for genetic background at centers

import os
import warnings
import pandas as pd

#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')
#%%

def genetic_background(data: str = "knockout"): 
    '''
    gets the genetic background of the controls and knockouts at each center
    
    import modules:
        import os
        import warnings
        import pandas as pd
    
    download gene and/or control data first
    
    Parameters
    ----------
     - data: str, uses either knockout or control data for strain access
         options: "knockout" (default), "control"
         optional prameter
    
    Returns
    -------
    None
    automatically assigns results data frame to variable
    '''
    
    # ignore warnings
    warnings.filterwarnings("ignore")
    
    ########
    #create folder if it doesn't already exist
    save = r'./overall_analysis'
    # check whether directory already exists, otherwise create new folder
    if not os.path.exists(save):
        os.mkdir(save)
    
    #create subfolder for saving if it doesn't already exist
    save = f'{save}/genetic_background'
    # check whether directory already exists, otherwise create new folder
    if not os.path.exists(save):
        os.mkdir(save)
    #########
    
    
    if data == "knockout":
        print("access genetic background with knockout data")
        read = r"./genes_all_columns/all_genes_all_columns.csv"
        df = pd.read_csv(read)
          
    else:
        print("access genetic background with control data")   
        read = r"./controls/recent_controls.csv"
        df = pd.read_csv(read)
    
    
    df = df.sort_values(by=["phenotyping_center"]).reset_index(drop=True)
    
    genetic_background = df.groupby(["phenotyping_center", "strain_name"]).size().reset_index(name='count')
        #the reset_index function is used to reset the index of the resulting DataFrame and assign a new 
        #column name "count" to the column that stores the group sizes ---> new 'count' column formed

    genetic_background.to_csv(f"{save}/genetic_background_{data}.csv", index=False)       
    globals()[f"genetic_background_{data}"] = genetic_background
    
#%%
genetic_background()
#gene or control data with all columns necessary