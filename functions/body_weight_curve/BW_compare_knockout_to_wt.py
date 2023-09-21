# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 17:34:20 2023

@author: Leonard Borchardt
"""


import os
import warnings
import pandas as pd
from datetime import datetime, timedelta

#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')
#%%

def BW_compare_knockout_to_wt(gene_symbols: list, parameter_stable_id: str = "IMPC_BWT_008_001", save_memory: bool = True):
    
    '''
    
    inserts in the gene_bw data frame in the cwd a new column with relative values for each body weight measurement to the controls of a time frame of +-14 days
    
    modules needed:
        import os
        import warnings
        import pandas as pd
        from datetime import datetime, timedelta
    
    
    Parameters
    ----------
    gene_symbol: list of strings, gene_symbols from IMPC
    parameter_stable_id: str, IMPC pipline parameter id that is connected to body weight curve obtaining, default: "IMPC_BWT_008_001"
    save_memory: bool, no returns, only saving, default: True, optional

    Returns (and saves):
    -------
    None
    automatically assigns variable for data frame with relative measurement values (relative to surrounding controls of the same sex)
    
    additional saving of the data frame
    '''
    
    # ignore warnings
    warnings.filterwarnings("ignore")
    
    print(parameter_stable_id)
    for gene_symbol in gene_symbols:
        print(gene_symbol)
        
        ########
        #create folder if it doesn't already exist
        path = f'./parameters/{parameter_stable_id}'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(path):
            os.mkdir(path)
        
        #create subfolder for saving if it doesn't already exist
        path = f'{path}/{gene_symbol}'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(path):
            os.mkdir(path)
        #########
        
        #read in and filter gene data    
        try:
            gene_df =  pd.read_csv(f"./genes/{gene_symbol}/{gene_symbol}.csv")
        except FileNotFoundError:
            print("download gene data first with function: download_genes.py")
            continue
        
        gene_df = gene_df[(gene_df.parameter_stable_id == parameter_stable_id)]
          
        #sort the data after specimen with ascending age
        gene_df = gene_df.sort_values(by = ["age_in_weeks"])
        gene_df = gene_df.sort_values(by = ["specimen_id"])
        gene_df = gene_df.reset_index(drop = True)
            #reset index without keeping old index column
        
        # make measured values positive
        gene_df['data_point'] = gene_df['data_point'].abs()
        
        
        #read in control data
        try:
            control_df = pd.read_csv(f"{path}/{gene_symbol}_{parameter_stable_id}_controls.csv")
        except FileNotFoundError:
            try:
                control_df = pd.read_csv(f"./parameters/{parameter_stable_id}/{parameter_stable_id}_controls.csv")
            except FileNotFoundError:
                print("Missing data frame with controls.")
                print("Run get_control_data function first!")
          
        # make measured values positive
        control_df['data_point'] = control_df['data_point'].abs()
          
        
        data_frame_list = []
        
        for specimen in set(gene_df.specimen_id):
            
            #phenotyping center of the animal
            center = gene_df[gene_df.specimen_id == specimen]["phenotyping_center"].reset_index(drop=True)[0]
            
            #get birthday of the animal
            date = gene_df[gene_df.specimen_id == specimen].date_of_birth.values
            date = date[0]
            
            #get surrounding controls
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
            dates = [date.strftime(date_format) for date in date_list]
            #that includes the birthday of the knockout +-14 surrounding days
            
            #filter control_df for controls in the time around knockout animal at the same phenotyping center
            control_data_dates = control_df[(control_df.phenotyping_center == center) & (control_df.date_of_birth.isin(dates))]
            
            ###### weekly measurements, unique for "body weight curve" parameter
            #make control data mean for every specimen age (in weeks)
            for week in set(gene_df[gene_df.specimen_id == specimen]["age_in_weeks"]):
                c_week_mean = control_data_dates[control_data_dates.age_in_weeks == week].data_point.mean()
                #control data mean for the particular age
                
                #reffering knockout measurement to that value
                rel_value = gene_df[(gene_df.specimen_id == specimen)&(gene_df.age_in_weeks == week)].data_point/c_week_mean
                
                d = pd.DataFrame({
                    "specimen_id": specimen,
                    "age_in_weeks": week,
                    "knockout_wt_ratio": rel_value
                    })
                
                data_frame_list.append(d)
    
        try:
            #concat the list of data frames with the specimen and week results
            rel_data = pd.concat(data_frame_list, ignore_index=True)
            
        except ValueError:
            print("missing knockout or control values for this parameter")
            continue
        
        #insert the rel. values in the original gene data frame for this parameter
        gene_df.insert(14, 'knockout_wt_ratio', rel_data['knockout_wt_ratio'])    
              #inserts it after equal index, so this has to rise continuously in both data frames for it to make sense  
        
        gene_df.to_csv(f'{path}/{gene_symbol}_{parameter_stable_id}_kwr.csv', index=False)
        
        if save_memory == False:
            globals()[f"{gene_symbol}_{parameter_stable_id}_kwr"] = gene_df
            #automatic variable assignment for the specific controls
            #no need for return statement
#%%

BW_compare_knockout_to_wt(["Prkab1", "Dbn1", "Ap4e1", "Nxn"])