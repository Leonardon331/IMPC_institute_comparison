# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 10:36:37 2023

@author: Leonard Borchardt
"""

## concatenate all knockout-control-ratios for all parameter_stable_ids

import os
import warnings
import pandas as pd

#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')

#%%


def concat_knockout_wt_ratios(BW_age_weeks: int = 14):
    '''
    import modules:
        import os
        import warnings
        import pandas as pd
    
    produces and saves concatenated kwr data frames under cwd/overall_analysis/knockout_wt_ratio
    produces a all ratios and a cleaned data frame with removed NaNs and double measurements and for IMPC_BWT_008_001 only the set age (see "Parameters")

    Parameters:
        BW_age_weeks: int, specify for Body weight curve (IMPC_BWT_008_001) series parameter the age in weeks, default: 14, optional
    -------
    Returns: None
    
    '''
    
    # ignore warnings
    warnings.filterwarnings("ignore")
    
    save = r'./overall_analysis'
    # check whether directory already exists, otherwise create new folder
    if not os.path.exists(save):
        os.mkdir(save)
    
    save = r'./overall_analysis/knockout_wt_ratio'
    # check whether directory already exists, otherwise create new folder
    if not os.path.exists(save):
        os.mkdir(save)
    
    contains = 'kwr'
    contains_not = "plot"         
    
    
    def add_kwrs(path = r"./parameters", kwrs = []):
        #recursive fct., append kwrs to a list
        for item in os.listdir(path):
            item_path = os.path.join(path, item)
            if os.path.isfile(item_path) and contains in item and contains_not not in item:
                try:
                    df = pd.read_csv(item_path)
                    kwrs.append(df)
                    print(item, "appended")
                except UnicodeDecodeError:
                    continue
                
            if os.path.isdir(item_path):
                add_kwrs(item_path, kwrs)
        return kwrs
    
    kwrs = add_kwrs()
    
    if len(kwrs) != 0:   
        #concatenate list to make data frame
        kwrs = pd.concat(kwrs,ignore_index=True)
        
        kwrs = kwrs.sort_values(by=["gene_symbol","parameter_stable_id","zygosity","sex", "phenotyping_center","specimen_id"]).reset_index(drop=True)
        kwrs.to_csv(f'{save}/all_knockout_wt_ratios.csv', index=False)
        print("all knockout-wt-ratios results concatenated")
        
        ### data cleaning
        #drop nan values
        kwrs_cleaned = kwrs.dropna(subset=['knockout_wt_ratio'])
        
        # incl. only week 14 for "IMPC_BWT_008_001" series parameter
        kwrs_cleaned.drop(kwrs_cleaned[(kwrs_cleaned.parameter_stable_id == "IMPC_BWT_008_001")&(kwrs_cleaned.age_in_weeks != BW_age_weeks)].index, inplace = True)
        
        ##remove multiple measurements
        ###from here: series parameters like bodyweight make not so much sense anymore### # -> set BW to set week (default: 14)
        #by forming the mean of data_point for each specimen (---> 1 data point for each specimen)
        kwrs_cleaned = kwrs_cleaned.groupby(["gene_symbol","procedure_stable_id","procedure_name","parameter_stable_id","parameter_name","specimen_id","phenotyping_center","sex","zygosity"])['knockout_wt_ratio'].mean().to_frame().reset_index()
            #to_frame() bcs. output is a Series, reset_index() to not have knockout_wt_ratio as index
        
        kwrs_cleaned = kwrs_cleaned.sort_values(by=["gene_symbol","parameter_stable_id","zygosity","sex", "phenotyping_center"]).reset_index(drop=True)
        kwrs_cleaned.to_csv(f'{save}/cleaned_knockout_wt_ratios.csv', index=False)
        print("results cleaned")
        
    else:
        print('check if the kwr data frames exist in parameters')
        
#%%
concat_knockout_wt_ratios()