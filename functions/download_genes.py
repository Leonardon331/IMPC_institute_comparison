# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 14:26:23 2023

@author: Leonard Borchardt
"""

import os
import warnings
import pandas as pd


#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')
#%%
def download_genes(gene_symbols: list):
    
    '''
    download and save gene data and save it as csv
    
    modules needed:
        import os
        import warnings
        import pandas as pd
        
    
    Parameters:
    gene_symbols: list of strings, gene_symbols for genes that should be downloaded and saved as csv files
    
    Returns:
        None
    '''
    
    # ignore warnings
    warnings.filterwarnings("ignore")
    
    for gene_symbol in gene_symbols:
        
        print(f"downloading {gene_symbol}")
        
        ########
        #create folder if it doesn't already exist
        path = r'./genes'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(path):
            os.mkdir(path)
        
        #create subfolder for saving if it doesn't already exist
        path = f'{path}/{gene_symbol}'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(path):
            os.mkdir(path)
        #########
        
        
        ########
        #create folder if it doesn't already exist
        path2 = r'./genes_all_columns'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(path2):
            os.mkdir(path2)
        #########
        
        #check if file already exists
        if os.path.exists(f'{path}/{gene_symbol}.csv'):
            print(f"{gene_symbol} already downloaded")
            continue
        
        
        url = f"https://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=gene_symbol:{gene_symbol}&wt=csv&rows=10000000"
        df = pd.read_csv(url)
        df = df.sort_values(by = ["procedure_name","parameter_stable_id","zygosity","sex","phenotyping_center"]).reset_index(drop=True)
        df.to_csv(f'{path2}/{gene_symbol}_all_columns.csv', index=False)
        print(f"{gene_symbol} downloaded")
        
        #url = f"https://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=gene_symbol:{gene_symbol}&fl=gene_symbol,genetic_background,procedure_name,procedure_stable_id,parameter_name,parameter_stable_id,zygosity,sex,phenotyping_center,biological_sample_group,colony_id,specimen_id,observation_type,data_point,category,metadata,date_of_experiment,date_of_birth,age_in_weeks,age_in_days,developmental_stage_name&wt=csv&rows=10000000"
        df = df[["gene_symbol","strain_name","procedure_name","procedure_stable_id","parameter_name","parameter_stable_id","zygosity","sex","phenotyping_center","biological_sample_group","colony_id","specimen_id","observation_type","data_point","category","metadata","date_of_experiment","date_of_birth","age_in_weeks","age_in_days","developmental_stage_name"]]
        df.to_csv(f'{path}/{gene_symbol}.csv', index=False)
        #globals()[f"{gene_symbol}"] = df

#%%
download_genes(["Prkab1", "Dbn1", "Ap4e1", "Nxn"])
