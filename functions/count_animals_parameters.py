# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 22:37:10 2023

@author: Leonard Borchardt
"""

import os
import warnings
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from num2words import num2words

#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')
#%%

def count_animals_parameters(gene_symbols: list):
    '''
    for given genes: counts the animals tested and parameters carried out at each phenotyping center

    import modules:
        import os
        import warnings
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns
        from num2words import num2words

    Parameters
    ----------
    gene_symbols: list of strings, symbols for genes

    Returns
    -------
    None
     - automatically assigns variables for gene data frames with phenotyping centers, animal parameter and animal count
     - saves the data framen in a certain directory in the cwd

    '''
    
    # ignore warnings
    warnings.filterwarnings("ignore")
    
    for gene_symbol in gene_symbols:
        
        print(gene_symbol)
        
        ########
        #create folder if it doesn't already exist
        path = r'./overall_analysis'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(path):
            os.mkdir(path)
        
        #create subfolder for saving if it doesn't already exist
        path = f'{path}/genes_counts'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(path):
            os.mkdir(path)
        
        #create subfolder for saving if it doesn't already exist
        path = f'{path}/{gene_symbol}'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(path):
            os.mkdir(path)
        #########
        
        #########
        
        #read in and filter gene data    
        df =  pd.read_csv(f"./genes/{gene_symbol}/{gene_symbol}.csv")
        
        
        parameters_centers = []
        # create list of parameters that were carried out at the center
        for center in set(df.phenotyping_center):       
            
            
                parameters = list(set(df[df.phenotyping_center == center]["parameter_stable_id"]))
                
                parameters = pd.DataFrame({
                    "phenotyping_center": center,
                    "parameter_stable_id": parameters
                    })
                
                parameters_centers.append(parameters)
        
        # create list of parameters that were carried out at the center
        parameters_centers = pd.concat(parameters_centers, ignore_index=True)
        parameters_centers = parameters_centers.sort_values(by=['phenotyping_center']).reset_index(drop=True)
          
        
        
        counts = []
        for center in set(df.phenotyping_center): 
            #number of animals and parameters for each center that tested the gene
            a_count = len(set(df[df.phenotyping_center == center]["specimen_id"]))
            p_count = len(set(df[df.phenotyping_center == center]["parameter_name"]))   ##parameter name bcs. of multiple ids for same parameter
            
            
            count = pd.DataFrame({
                "phenotyping_center": center,
                "animal_count": [a_count],
                "parameter_count": [p_count]
                })
            
            counts.append(count)
            
        
        #number of animals and parameter overlap across all centers
        a_count = [len(set(df.specimen_id)), float("nan"), float("nan")]
            #all animals used for that gene
            
        
        #parameter overlap across centers (--> mandatory IMPC measurements!, but prob not all...)
        all_phenotyping_centers = pd.Series(list(set(df.phenotyping_center)))
        center_count = len(all_phenotyping_centers)
        
        print("check parameter overlap between all centers")
        cross_center_ps_ids = []
        most_center_ps_ids = []
        df_gr = df.groupby("parameter_stable_id")
        for name,group in df_gr:
            print(f'checking {gene_symbol}, {name} for center overlap')
            s = all_phenotyping_centers.isin(set(group.phenotyping_center)) #Boolean Seeries
            id_num_centers = len(s.loc[lambda x : (x == True)])
            print(f'{name} at {id_num_centers} center(s)')
            if all(s): # if all phhenotyping centers in ps_id group phenotyping centers True append ps_id
                cross_center_ps_ids.append(name)
                print(f"{name} measured at every center!")
            if id_num_centers >=center_count-2:
                most_center_ps_ids.append(name)
                    
                    
        p_count = [len(set(df.parameter_name)),len(cross_center_ps_ids),len(most_center_ps_ids)] #count center overarching parameters and most one, too
        

        count = pd.DataFrame({
            "phenotyping_center": "all",
            "animal_count": [a_count[0]],
            "parameter_count": [p_count[0]]
            })
        
        counts.append(count)
        
        count = pd.DataFrame({
            "phenotyping_center": "overlap-all",
            "animal_count": [a_count[1]],
            "parameter_count": [p_count[1]]
            })
        
        counts.append(count)

        
        count = pd.DataFrame({
            "phenotyping_center": f"overlap-{num2words(center_count-2)}+",
            "animal_count": [a_count[2]],
            "parameter_count": [p_count[2]]
            })
        
        counts.append(count)

        

        #### save also id lists
       
        cross_center_ps_ids = pd.DataFrame(dict(
            parameter_stable_id = cross_center_ps_ids
            ))
        
        cross_center_ps_ids = cross_center_ps_ids.sort_values(by=['parameter_stable_id']).reset_index(drop=True)
        cross_center_ps_ids.to_csv(f"./genes/{gene_symbol}/{gene_symbol}_cross_center_ps_ids.csv", index = False)
        globals()[f"{gene_symbol}_cross_center_ps_ids"] = cross_center_ps_ids
        
       
        most_center_ps_ids = pd.DataFrame(dict(
            parameter_stable_id = most_center_ps_ids
            ))
    
        most_center_ps_ids = most_center_ps_ids.sort_values(by=['parameter_stable_id']).reset_index(drop=True)
        most_center_ps_ids.to_csv(f"./genes/{gene_symbol}/{gene_symbol}_most_center_ps_ids.csv", index = False)
        globals()[f"{gene_symbol}_most_center_ps_ids"] = most_center_ps_ids
        
        ####
        
        
        #concatenate counts list of count data frames
        counts = pd.concat(counts, ignore_index=True)
        counts = counts.sort_values(by=['phenotyping_center']).reset_index(drop=True)
            #How does sorting the function work?
        counts.to_csv(f'{path}/{gene_symbol}_counts.csv',index=False)
        globals()[f"{gene_symbol}_counts"] = counts

        #####make bar plot
        # Create subplots with adjusted figsize
        #adjust the subplot layout
        subplot_height = 8  # adjust this value to control the height
        subplot_width = 8   # adjust this value to control the width
        
        subplot_rows = 1
        subplot_columns = 2
        
        total_fig_height = subplot_height * subplot_rows
        total_figure_width = subplot_width * subplot_columns
        
        # Create subplots with adjusted figsize
        fig, axs = plt.subplots(subplot_rows, subplot_columns, figsize=(total_figure_width, total_fig_height), dpi=500)
    
        sns.set_theme()
                                    #slicing stop not incl.
        sns.barplot(data=counts[0:len(counts)-3], x="phenotyping_center", y="animal_count", ax=axs[0])    #data=counts[0:len(counts)-2] to exclude categories with animal_count nan values + all animals
        axs[0].set_ylabel("Animal Count", fontsize=14)
        
        
        # Change the label "all" to "all-distinct" in the counts DataFrame
        counts.loc[counts['phenotyping_center'] == 'all', 'phenotyping_center'] = 'all-distinct'
            # df.loc[row, column_name]
        
        sns.barplot(data=counts[counts.phenotyping_center!="all-distinct"], x="phenotyping_center", y="parameter_count", ax=axs[1])
        axs[1].set_ylabel("Parameter Count", fontsize=14)
        
        fig.subplots_adjust(top=.9)
        fig.suptitle(f'Number of animals and parameters across centers for {gene_symbol}', size = 18, ha = 'center')
        #bug: ha/horizontalaligment= 'right' ----> 'left'
            
        for ax in axs.flat:
            ax.tick_params(axis='x', labelrotation = 90)
            ax.set_xlabel("Phenotyping Center", fontsize=14)
            
            for i in ax.containers:
                ax.bar_label(i,)
        
        plt.savefig(f"{path}/bar_plot_{gene_symbol}_counts.png", dpi = 500, bbox_inches='tight')
        #whole path for plot saving important
        
        ##############
           
        ###
        #list of all parameters measured
        parameters_centers.to_csv(f'./genes/{gene_symbol}/{gene_symbol}_all_parameters_at_centers.csv',index=False)
        ###
    
        print("--------------------------------------------------")
        print()
#%%
#call the function

count_animals_parameters(gene_symbols = ["Prkab1", "Dbn1", "Ap4e1", "Nxn"])