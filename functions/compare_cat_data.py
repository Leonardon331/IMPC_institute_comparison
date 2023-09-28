# -*- coding: utf-8 -*-
"""
Created on Sun Sep  3 10:13:56 2023

@author: Leonard Borchardt
"""

import os
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')
#%%
def compare_cat_data(gene_symbols: list, min_centers: int = 4, focus_zygosity: bool = True, skip_plotting: bool = False, plot_only_differences: bool = True):
    '''
    
    plots for categorical data the proportion between the categories in a bar plot
    for each parameter (that has much data across centers) with bars symbolizing the centers
       - saves the plot under cwd/overall_analysis/cat_data/"gene_symbol"
       - saves data frame with category proportions in the same directory
    
    run download_genes.py first
    
    import modules:
        import os
        import warnings
        import pandas as pd
        import numpy as np
        import matplotlib.pyplot as plt
        import seaborn as sns
    
    Parameters:
        gene_symbols: list of strings, gene symbols for genes downloaded
        min_centers: int, min. phenotyping centers for each comparison, default: 4, increasing leads to less memory consumption, optional
        focus_zygosity: bool, only zygosity with the most data across phenotyping centers is used for comparision, default: True, optional
        skip_plotting: bool, skip plotting and only calculate category proportions, saves memory, defaut: False, optional
        plot_only_differences: bool, plot only interesting differences between the institutes, default: True, optional
        
    Returns: None
    
    Saves: data frame with category proportions for each center and parameter
    '''
    
    # ignore warnings
    warnings.filterwarnings("ignore")
    
    print()
    print()
    
    for gene_symbol in gene_symbols:
        print(gene_symbol)
        read = f'./genes_all_columns/{gene_symbol}_all_columns.csv'   ##read in the gene data
        df = pd.read_csv(read)
    
        ########
        #create folder if it doesn't already exist
        save = r'./overall_analysis'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(save):
            os.mkdir(save)
        
        #create subfolder if it doesn't already exist
        save = f'{save}/categorical_data'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(save):
            os.mkdir(save)
        
        #create subfolder if it doesn't already exist
        save = f'{save}/{gene_symbol}'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(save):
            os.mkdir(save)
            
        #########
        
        # only cat. data
        df = df[df.observation_type == "categorical"]
        
        #drop nan values
        df = df.dropna(subset=['category'])
        
        
        #### additional filtering if parameter focus_zygosity = True
        if focus_zygosity == True:
            #choose biggest zygosity subdataset for center comparison
            if df[df.zygosity == "homozygote"].size > df[df.zygosity == "heterozygote"].size and df[df.zygosity == "homozygote"].size > df[df.zygosity == "hemizygote"].size:
                zygosity = "homozygote"
                df = df[df.zygosity == zygosity]
                
            elif df[df.zygosity == "heterozygote"].size > df[df.zygosity == "hemizygote"].size:
                zygosity = "heterozygote"
                df = df[df.zygosity == zygosity]
            
            else:
                zygosity = "hemizygote"
                df = df[df.zygosity == zygosity]
          ####
          
        
        #relative amount of categories
        rel_categories =  []
        
        parameters = df.groupby("parameter_stable_id")
        for name, parameter in parameters:
            
            center_zyg_sex = parameter.groupby(["phenotyping_center","zygosity", "sex"])
            for name_, group in center_zyg_sex:
                all_cat_count = group.size
                for category in set(group.category):
                    cat_count = group[group.category == category].size
                    rel_cat_count = cat_count/all_cat_count
                    
                    rel_category = pd.DataFrame(dict(
                        gene_symbol = [gene_symbol],
                        parameter_stable_id = [name],
                        parameter_name = [parameter.parameter_name.reset_index(drop=True)[0]],
                        phenotyping_center = [group.phenotyping_center.reset_index(drop=True)[0]],
                        zygosity = [group.zygosity.reset_index(drop=True)[0]],
                        sex = [group.sex.reset_index(drop=True)[0]],
                        category = [category],
                        cat_count = [cat_count],
                        cat_proportion = [rel_cat_count]
                        ))
                    
                    rel_categories.append(rel_category)
                    
                    
        rel_categories = pd.concat(rel_categories, ignore_index=True)
        
        del df
        # save memory
        
        #remove centers with less than 3 samples for meaningful statistics
        groups = []
        rel_categories_groups = rel_categories.groupby(["parameter_stable_id", "phenotyping_center", "zygosity", "sex", "category"])
        for name,group in rel_categories_groups:
            if group.cat_count.values[0] < 3:
                continue
            else:
                groups.append(group)       
        if len(groups) != 0:
            rel_categories = pd.concat(groups,ignore_index=True)
        else:
            print(f"not enough data for {gene_symbol} for categorical center comparison")
            print()
            continue    #next gene
        
        #save data frame
        rel_categories = rel_categories.sort_values(by=['parameter_stable_id', 'zygosity', 'sex', 'phenotyping_center']).reset_index(drop=True)
        rel_categories.to_csv(f'{save}/{gene_symbol}_category_proportions.csv', index=False)
        print("relative center values for categories calculated")
        
        if skip_plotting == True:
            print("------------------------")
            continue
        
        
        
        #### make bar plots ###
        print("make bar plots for each parameter")
        
        ## for stacked par plots it's important to have the proportions in different columns
        ## ---> create for every parameter a new data frame with categories as columns
        
        #remove groups with less than min_centers for meaningful plots
        groups = []
        rel_categories_groups = rel_categories.groupby(["parameter_stable_id", "sex", "zygosity"])
        for name,group in rel_categories_groups:
            if len(set(group.phenotyping_center)) < min_centers:
                continue
            else:
                groups.append(group)  
        if len(groups) != 0:
            rel_categories = pd.concat(groups,ignore_index=True)
        else:
            print(f"not enough data for {gene_symbol} not enough for categorical center comparison")
            print()
            continue    #next gene

        
        # plots for each parameter separately
        parameters = list(set(rel_categories.parameter_stable_id))
        parameters.sort()
        
        
        for parameter in parameters:
            print(parameter)
            if os.path.exists(f'{save}/barplot_{gene_symbol}_{zygosity}_{parameter}_rel_categories.png'):
                print("bar plot already created")
                continue
            
            data = rel_categories[rel_categories.parameter_stable_id == parameter].reset_index(drop=True)
                
            categories = data.sort_values(by = ["cat_count"]).category.unique().tolist()
            #print(categories)
            
            if plot_only_differences == True:
                #plot only results, where there is a interesting difference between the centers
                for category in categories:
                    if all((data.cat_proportion == 1) & (data.category == category)):
                        x = 0
                        break   
                    else:
                        x = 1
                if x == 0:
                    continue   #next parameter
            
            sexes = data.sort_values(by = ["sex"], ascending = False).sex.unique().tolist()
            
            parameter_name = data.parameter_name[0]
            
            zygosities = list(set(data.zygosity))
            zygosities.sort()
            
            centers = list(set(data.phenotyping_center))
            centers.sort()
            
            d = pd.DataFrame(columns=["phenotyping_center","zygosity","sex"] + categories)
            for center in centers:
                for zygosity in zygosities:
                    for sex in sexes:
                        row_list = [center, zygosity, sex]      # to be filled and later appended to d
                        d_row = data[(data.phenotyping_center == center)&(data.zygosity == zygosity)&(data.sex == sex)]
                        for category in categories:
                            if category in list(d_row.category):
                            #in list(Series) or Series.isin()
                                row_list.append(d_row[d_row.category == category].cat_proportion.reset_index(drop=True)[0].item())
                                #append category proportion of this parameter/center/zygosity/sex
                            else:
                                row_list.append(0)  #category for this parameter/center/zygosity/sex not observed
                                
                        #append row_list to d       
                        d.loc[len(d)] = row_list
            
            
            ##reorder columns
            # Calculate column sums and reorder columns in descending order
            column_sums = d.select_dtypes(include=[np.number]).sum().sort_values(ascending=False)
            # Create a list of column names sorted by descending sum
            sorted_columns = column_sums.index.tolist()
                #.index to get column names
            #sorted column names:
            sorted_columns = ["phenotyping_center","zygosity","sex"] + sorted_columns
            # reorder DataFrame columns based on sorted column names
            d = d[sorted_columns]
            
            #adjust the subplot layout
            subplot_height = 8  # adjust this value to control the height
            subplot_width = 8   # adjust this value to control the width
            
            subplot_rows = len(zygosities)
            subplot_columns = len(sexes)
            
            total_fig_height = subplot_height * subplot_rows
            total_figure_width = subplot_width * subplot_columns
            
            # Create subplots with adjusted figsize
            fig, axs = plt.subplots(subplot_rows, subplot_columns, figsize=(total_figure_width, total_fig_height), dpi=500)
           
            sns.set_theme()              
            
            if len(sexes) == 1:
                #no iteration over axs necessary/working, only one axis
                d.plot(kind='bar', stacked=True, ax = axs)
                axs.set_xticklabels(d.phenotyping_center)
                axs.tick_params(axis='x', labelrotation = 90)
                axs.set_title(sexes[0].capitalize(), fontsize = 14)
                axs.set_xlabel("Phenotyping Center", fontsize=14)
                axs.set_ylabel("Proportion", fontsize=14)
                axs.get_legend().remove()  
                
            else:
                i = 0
                for ax in axs.flat:
                    d_sex=d[d.sex == sexes[i]]
                    d_sex.plot(kind='bar', stacked=True, ax=ax)
                    ax.set_xticklabels(d_sex.phenotyping_center)
                    ax.tick_params(axis='x', labelrotation = 90)
                    ax.set_title(sexes[i].capitalize(), fontsize = 14)
                    ax.set_xlabel("Phenotyping Center", fontsize=14)
                    ax.set_ylabel("Proportion", fontsize=14)
                    ax.get_legend().remove()
                    i+=1

            #legend
            plt.legend(bbox_to_anchor=(1.02,1), loc='upper left', title = "Category")
            
            
            fig.subplots_adjust(top=.85)
            fig.suptitle(f'{gene_symbol}/{zygosity}: Proportion of Categories for {parameter_name} ({parameter})', size = 18, ha = 'center')
            plt.tight_layout()
            
            plt.savefig(f'{save}/barplot_{gene_symbol}_{zygosity}_{parameter}_rel_categories.png', dpi = 500, bbox_inches='tight')

            plt.clf() #clear current figure, save memory
            
        del rel_categories #free memory  
        print("------------------------")
#%%
gene_symbols = ["Ap4e1","Dbn1","Nxn","Prkab1"]
compare_cat_data(gene_symbols, min_centers = 1, skip_plotting=True)