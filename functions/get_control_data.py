# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 11:27:27 2023

@author: Leonard Borchardt
"""

import os
import pandas as pd
from urllib.error import HTTPError, URLError
from datetime import datetime, timedelta

#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')
#%%
def get_control_data(gene_symbols: list, parameter_stable_ids: list, save_memory: bool = True):
    
    '''
    returns a data frame for the controls surrounding the birthday of particular knockout mice used in a certain parameter test across centers for +-14 days
     - if no timewise close control data available, try to get time independent control data of each center
    
    Internet connection needed!
    csv-gene-dataframe in cwd needed (gene_symbol.csv)
    ## running time drastically with a lot with longer function parameter lists
    
    modules needed:
        import os
        import pandas as pd
        from urllib.error import HTTPError, URLError
        from datetime import datetime, timedelta

    Parameters
    ----------
    gene_symbols: list of strings, gene symbols for knocked-out gene
    parameter_stable_ids: list of strings, phenotyping parameter ids for different parameter measurements
    save_memory: bool, no assigned variables only saving of the data frame

    Returns
    -------
    None
    automatic variable assignment: 
        control mice dataframe that where phenotyped for this parameter timewise close to knockout mice
    '''
    
    for gene_symbol in gene_symbols:
        
        #read in gene_data from:
        try:
            df = pd.read_csv(f"genes/{gene_symbol}/{gene_symbol}.csv")
        except FileNotFoundError:
            print("download gene data first with function: download_genes.py")
            continue
        
        for parameter_stable_id in parameter_stable_ids:
            print(gene_symbol)
            print(parameter_stable_id)

            file = f'./parameters/{parameter_stable_id}/{gene_symbol}/{gene_symbol}_{parameter_stable_id}_controls.csv'
            if os.path.exists(file):
                print("already downloaded")
                print("--------------------")
                continue
            
            ######## write to:
            #create folder if it doesn't already exist
            path = f'./parameters/{parameter_stable_id}'
            # check whether directory already exists, otherwise create new folder
            if not os.path.exists(path):
                os.mkdir(path)
            
            #create subfolder for saving if it doesn't already exist
            path = f'{path}/{gene_symbol}'
            # check whether directory already exists (if yes, then move to next iteration), otherwise create new folder
            if not os.path.exists(path):
                os.mkdir(path)    
            #########
            
            
            #create empty WT data frame
            WT = pd.DataFrame()
            
            #create empty list to append the wt data frames later
            wt_list = []
            
            #get WT data for same phenotyping_center near the knockout specimen birthday (--> "hard window")
            for center in set(df["phenotyping_center"]):
                
                #birthdays of knockout mice at the center
                dates = set(df[df.phenotyping_center == center]["date_of_birth"].dropna())
                dates = list(dates)
                    #conversion to list to be able to append values later
                #----------------------
                #Add more dates to access more of control data for more reliable baseline
                #add dates 14 days in one direcion and 14 days in the other direction from knockout specimen birthday
                
                added_days = []
                
                for date in dates:
                    
                    # Convert the initial date string to a datetime object
                    date_format = '%Y-%m-%dT%H:%M:%SZ'
                    current_date = datetime.strptime(date, date_format)
                    
                    # Create an empty list to store the dates
                    date_list = [current_date]
                    
                    # Define the number of dates you want to add
                    num_dates_to_add = 14
                    
                    # Append the following dates to the list, + 14 days
                    for _ in range(num_dates_to_add):
                        # Add one day to the current date
                        current_date = current_date + timedelta(days=1)
                        date_list.append(current_date)
                    
                    current_date = datetime.strptime(date, date_format)
                    # Append the following dates to the list, - 14 days
                    for _ in range(num_dates_to_add):
                        # Add one day to the current date
                        current_date = current_date + timedelta(days=-1)
                        date_list.append(current_date)
                    
                    # Convert the list of datetime objects back to formatted strings
                    formatted_date_list = [date.strftime(date_format) for date in date_list]
                    
                    for day in formatted_date_list:
                        added_days.append(day)
                
                for day in added_days:
                    dates.append(day)
                
                dates.sort()
                #unique dates
                dates = set(dates)
                #now for every knockout birthday of the the center, 14 others in both directions (---> more controls accessable)
                
                #--------------------------
                #get WT data
            
                for date in dates:
                    try:
                        url = f'https://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=biological_sample_group:control%20AND%20parameter_stable_id:{parameter_stable_id}%20AND%20phenotyping_center:"{center}"%20AND%20date_of_birth:"{date}"&fl=biological_sample_group,strain_name,phenotyping_center,procedure_name,procedure_stable_id,parameter_name,parameter_stable_id,sex,colony_id,specimen_id,observation_type,data_point,category,metadata,date_of_experiment,date_of_birth,age_in_weeks,age_in_days,developmental_stage_namee&wt=csv&rows=10000000'
                        #get only control data for given parameter_stable_id (download takes otherwise too much time)    
                        #url needs center with spaces (UC Davis) in double quotation marks "", date needs to be surrounded by "" too
                        if " " in url:
                            print("replace space in url")
                            url = url.replace(" ", "%20")
            
                        wt = pd.read_csv(url)
                        
                        if wt.empty:
                            print(f"no control data for {gene_symbol} for {parameter_stable_id} at {center} at day {date}")
                            
                        else:
                            wt_list.append(wt)
                            print(f"WT data for {gene_symbol} for {parameter_stable_id} at {center} on {date} successfully collected!")
                            
                    except (HTTPError, URLError):  #TimeoutError, 
                    #check for multiple errors with tuple
                        print("url error, check maybe parameters")
                    
                        
            
            #concat list of data frames (wt_list) to one data frame (WT)
            try:
                WT = pd.concat(wt_list, ignore_index=True)
                WT = WT.sort_values(by=["procedure_name","parameter_stable_id", "sex","phenotyping_center"]).reset_index(drop=True)
                WT.to_csv(f'{path}/{gene_symbol}_{parameter_stable_id}_controls.csv', index=False)
                
                if save_memory == False:
                    globals()[f"{gene_symbol}_{parameter_stable_id}_controls"] = WT
                    #automatic variable assignment for the specific controls
                
                
            except ValueError:
                #if wt_list is still emty (pd.concat will raise ValueError)
                print("no closely measured control data")
                print("controls maybe not measured permanently for this parameter")
                print()
                print("get time independend control data of each center instead")
                
                for center in set(df["phenotyping_center"]):
                    try:
                        url = f'https://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=biological_sample_group:control%20AND%20parameter_stable_id:{parameter_stable_id}%20AND%20phenotyping_center:"{center}"&fl=biological_sample_group,strain_name,procedure_name,procedure_stable_id,parameter_name,parameter_stable_id,sex,phenotyping_center,colony_id,specimen_id,observation_type,data_point,category,metadata,date_of_experiment,date_of_birth,age_in_weeks,age_in_days,developmental_stage_name&wt=csv&rows=10000000'
                        ## url without day_of_birth query this time
                        #get only control data for given parameter_stable_id (download takes otherwise too much time)    
                        #url needs center with spaces (UC Davis) in double quotation marks "", date needs to be surrounded by "" too
                        if " " in url:
                            print("replace space in url")
                            url = url.replace(" ", "%20")
            
                        wt = pd.read_csv(url)
                        
                        if wt.empty:
                            print(f"also no time independent control data for {gene_symbol} for {parameter_stable_id} at {center}")
                           
                        else:
                            wt_list.append(wt)
                            print(f"Time independent WT data for {gene_symbol} for {parameter_stable_id} at {center} successfully collected!")
                            
                    except (HTTPError, URLError):   #(+TimeoutError)
                    #check for multiple errors with tuple
                        print("url error, check maybe parameters")
                    
                    
                #again try concatenation of list of data frames (wt_list) to one data frame (WT)
                try:
                    WT = pd.concat(wt_list, ignore_index=True)
                    WT = WT.sort_values(by=["procedure_name","parameter_stable_id", "sex","phenotyping_center"]).reset_index(drop=True)
                    WT.to_csv(f'./parameters/{parameter_stable_id}/{parameter_stable_id}_controls.csv', index=False)
                    
                    if save_memory == False:
                        globals()[f"{parameter_stable_id}_controls"] = WT
                        #automatic variable assignment for the specific controls
                    
                except ValueError:
                    print(f"no data for {gene_symbol} under {parameter_stable_id}")
                
            print("------------------------------------------------------")
            
            
#%%

#call function
gene_list = ["Prkab1", "Dbn1", "Ap4e1", "Nxn"]
parameters = ["IMPC_ABR_004_001",
"IMPC_ABR_006_001",
"IMPC_ABR_008_001",
"IMPC_ABR_010_001",
"IMPC_ABR_012_001",
"IMPC_ACS_001_001",
"IMPC_ACS_002_001",
"IMPC_ACS_003_001",
"IMPC_ACS_004_001",
"IMPC_ACS_006_001",
"IMPC_ACS_007_001",
"IMPC_ACS_008_001",
"IMPC_ACS_009_001",
"IMPC_ACS_033_001",
"IMPC_ACS_034_001",
"IMPC_ACS_035_001",
"IMPC_ACS_037_001",
"IMPC_CAL_001_001",
"IMPC_CAL_002_001",
"IMPC_CAL_017_001",
"IMPC_CBC_004_001",
"IMPC_CBC_005_001",
"IMPC_CBC_006_001",
"IMPC_CBC_007_001",
"IMPC_CBC_008_001",
"IMPC_CBC_009_001",
"IMPC_CBC_010_001",
"IMPC_CBC_012_001",
"IMPC_CBC_013_001",
"IMPC_CBC_014_001",
"IMPC_CBC_015_001",
"IMPC_CBC_016_001",
"IMPC_CBC_017_001",
"IMPC_CBC_018_001",
"IMPC_CSD_032_001",
"IMPC_DXA_001_001",
"IMPC_DXA_002_001",
"IMPC_DXA_003_001",
"IMPC_DXA_004_001",
"IMPC_DXA_005_001",
"IMPC_DXA_007_001",
"IMPC_DXA_008_001",
"IMPC_DXA_009_001",
"IMPC_DXA_010_001",
"IMPC_GRS_003_001",
"IMPC_GRS_008_001",
"IMPC_GRS_009_001",
"IMPC_GRS_010_001",
"IMPC_GRS_011_001",
"IMPC_HEM_001_001",
"IMPC_HEM_002_001",
"IMPC_HEM_003_001",
"IMPC_HEM_004_001",
"IMPC_HEM_005_001",
"IMPC_HEM_006_001",
"IMPC_HEM_007_001",
"IMPC_HEM_008_001",
"IMPC_HWT_007_001",
"IMPC_HWT_008_001",
"IMPC_HWT_012_001",
"IMPC_IPG_001_001",
"IMPC_IPG_010_001",
"IMPC_IPG_011_001",
"IMPC_IPG_012_001",
"IMPC_OFD_009_001",
"IMPC_OFD_010_001",
"IMPC_OFD_012_001",
"IMPC_OFD_014_001",
"IMPC_OFD_016_001",
"IMPC_OFD_020_001",
"IMPC_PAT_049_002"]
get_control_data(gene_list,parameters)