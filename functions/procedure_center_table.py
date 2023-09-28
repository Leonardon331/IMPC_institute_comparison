# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 14:03:14 2023

@author: Leonard Borchardt
"""


import os
import warnings
import pandas as pd
import re

#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')
#%%
def procedure_center_table(data: str = "knockout", procedures: list = ["VIA","FER","BWT","OFD","GRS","CSD","ACS","ECG","IPG","ABR","XRY","DXA","EYE","CBC","PAT","HEM","HWT","IMM"]):
    '''
    makes a table procedure id (procedure_stable_id) comparison across phenotyping centers for mandatory procedures
    
    import modules:
        import os
        import pandas as pd

    needs following data frame:
        for data == "knockout":
            ./genes_all_columns/all_genes_all_columns.csv
        for data == "control":
            ./controls/recent_controls.csv
        
        if not existing:
            - run for "knockout": download_genes.py, or for "control": get_recent_controls.py
            - run for "knockout" additionally: concat_all_genes.py

    Parameters
    ----------
    procedures : list of strings, default: mandatory IMPC procedures, optional
    data: str, data frame used for center-procedure table creation, values: "knockout" (default) or "control", optional

    Returns
    -------
    None.
    saves table in ./overall_analysis/procedures
    automatically assigns variable for the resulting data frame

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
    save = f'{save}/procedures'
    # check whether directory already exists, otherwise create new folder
    if not os.path.exists(save):
        os.mkdir(save)
    #########
    
    table = pd.DataFrame(columns = ["phenotyping_center"]+procedures)
    
    if data == "knockout":
        print("using knockout data for center-procedure table")
        df = pd.read_csv(r'./genes_all_columns/all_genes_all_columns.csv')
        
        df_gr = df.groupby("phenotyping_center")
        
        for name, group in df_gr:
            center_val = []
            center_val.append(name)
            for procedure in procedures:
                if group.procedure_stable_id.str.contains(procedure, regex = False).any():
                    #regex=True, would prob. also work
                    procedure_stable_id = group[group.procedure_stable_id.str.contains(procedure, regex = False)].procedure_stable_id.reset_index(drop=True)[0]
                    center_val.append(procedure_stable_id)
                    
                else:
                    center_val.append(float("nan"))
                
            #add "center_val" list as row to data frame "table"
            table.loc[len(table)] = center_val
        
        
        #transpose the data frame (while removing procedure index and changing the column names)
        table = table.transpose().reset_index(drop = True)
        table.columns = table.iloc[0]
        table = table.drop(table.index[0]).reset_index(drop = True)
        
        table.insert(0,"Procedure", procedures) 
    
    else:
        print("using recent control data for center-procedure table")
        df = pd.read_csv(r"./controls/recent_controls.csv")
        
        df_gr = df.groupby("phenotyping_center")
        
        for name, group in df_gr:
            center_val = []
            center_val.append(name)
            for procedure in procedures:
                #control data frame has only parameter_stable_id ---> create procedure_stable_id
                if group.parameter_stable_id.str.contains(procedure, regex = False).any():
                    #regex=True, would prob. also work
                    parameter_stable_id = group[group.parameter_stable_id.str.contains(procedure, regex = False)].parameter_stable_id.reset_index(drop=True)[0]
                    
                    pattern_to_remove = r'\d+_'  # This pattern matches one or more digits followed by an underscore

                    procedure_stable_id = re.sub(pattern_to_remove, '', parameter_stable_id)
                    #print(procedure_stable_id)
                    
                    center_val.append(procedure_stable_id)
                    
                else:
                    center_val.append(float("nan"))
                
            #add "center_val" list as row to data frame "table"       
            table.loc[len(table)] = center_val
        
        
        #transpose the data frame (while removing procedure index and changing the column names)
        table = table.transpose().reset_index(drop = True)
        table.columns = table.iloc[0]
        table = table.drop(table.index[0]).reset_index(drop = True)
        
        table.insert(0,"Procedure", procedures)
    
    table.to_csv(f"{save}/table_center_procedure_{data}.csv", index=False)
    globals()[f"table_center_procedure_{data}"] = table
    print(f"center-procedure table created with {data} data")

    
#%%
procedure_center_table()

