# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 23:14:36 2023

@author: Leonard Borchardt
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Sep  3 10:13:56 2023

@author: Leonard Borchardt
"""

#first plt bug, redo

import os
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')
#%%
def plot_cat_combined(gene_symbols: list, parameter_stable_ids: list, sexes: list, save_memory: bool = False):
    '''
    
    plots for categorical data the proportion between categories in a combined bar plot for all genes
    
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
        parameter_stable_ids: list of strings, IMPC parameters that should be plotted
        sexes: list of strings, "male" and "female" possible
        save_memory: bool, saves memory by not showing the plots (only saving), default: False, optional
        
    Returns: None
    
    Saves: combined plot with category proportions for each center and gene
    '''
    
    # ignore warnings
    warnings.filterwarnings("ignore")
    
    print()
    print()
    
    
    for parameter_stable_id in parameter_stable_ids:
        print(parameter_stable_id)
        for sex in sexes:
            print(sex)
            d_list = [] #append data frames in wide format to, later concatenation
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
                save = f'{save}/genes_combined'
                # check whether directory already exists, otherwise create new folder
                if not os.path.exists(save):
                    os.mkdir(save)
                    
                #create subfolder if it doesn't already exist
                save = f'{save}/{parameter_stable_id}'
                # check whether directory already exists, otherwise create new folder
                if not os.path.exists(save):
                    os.mkdir(save)
                    
                #########
                
                # only cat. data
                df = df[(df.observation_type == "categorical")&(df.parameter_stable_id==parameter_stable_id)&(df.sex==sex)]
                
                #drop nan values
                df = df.dropna(subset=['category'])
                
                
                #### additional filtering 
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
                
                #save data frame
                rel_categories = rel_categories.sort_values(by=['parameter_stable_id', 'zygosity', 'sex', 'phenotyping_center']).reset_index(drop=True)
                print("relative center values for categories calculated")
                
                
                
                #### make bar plots ###
                print("make bar subplot")
                
                ## for stacked par plots it's important to have the proportions in different columns
                ## ---> create for every parameter a new data frame with categories as columns
                    
                data = rel_categories
                    
                categories = data.category.unique().tolist()
                #print(categories)
                
                sexes = data.sort_values(by = ["sex"], ascending = False).sex.unique().tolist()
                
                parameter_name = data.parameter_name[0]
                
                zygosities = list(set(data.zygosity))
                zygosities.sort()
                
                centers = list(set(data.phenotyping_center))
                centers.sort()
                
                #bulild up data frame d with wide format
                d = pd.DataFrame(columns=["gene_symbol","phenotyping_center","zygosity","sex"] + categories)
                for center in centers:
                    for zygosity in zygosities:
                        for sex in sexes:
                            row_list = [gene_symbol, center, zygosity, sex]      # to be filled and later appended to d
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
                d_list.append(d)
                    
            d_concat = pd.concat(d_list, ignore_index=True)
                    
             
            ##reorder columns
            # Calculate column sums and reorder columns in descending order
            column_sums = d_concat.select_dtypes(include=[np.number]).sum().sort_values(ascending=False)
            # Create a list of column names sorted by descending sum
            sorted_columns = column_sums.index.tolist()
                #.index to get column names
            #sorted column names:
            sorted_columns = ["gene_symbol","phenotyping_center","zygosity","sex"] + sorted_columns
            # reorder DataFrame columns based on sorted column names
            d_concat = d_concat[sorted_columns]
             
            
            #adjust the subplot layout
            subplot_height = 8  # adjust this value to control the height
            subplot_width = 8   # adjust this value to control the width
            
            subplot_rows = 2
            subplot_columns = 2
            
            total_fig_height = subplot_height * subplot_rows
            total_figure_width = subplot_width * subplot_columns
            
            # Create subplots with adjusted figsize
            fig, axs = plt.subplots(subplot_rows, subplot_columns, figsize=(total_figure_width, total_fig_height), dpi=500)
               
            sns.set_theme()
            
            ax_list = []
            for ax in axs.flat:
                ax_list.append(ax)
               
            i = 0      
            for gene_symbol in gene_symbols: 
                d = d_concat[d_concat.gene_symbol==gene_symbol]
                zygosity = d.zygosity.to_numpy()[0]
                
                ax = ax_list[i]
                
                d.plot(kind='bar', stacked=True, ax=ax, legend = False)
                ax.set_xticklabels(d.phenotyping_center)
                ax.tick_params(axis='x', labelrotation = 90)
                ax.set_title(f"{gene_symbol}, {zygosity} {sex}s", fontsize = 18)
                ax.set_xlabel("Phenotyping Center", fontsize=14)
                ax.set_ylabel("Proportion", fontsize=14)
                i+=1
            
            #legend
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            fig.legend(handles=by_label.values(), labels=by_label.keys(), bbox_to_anchor=(1,0.95), loc='upper left', title = "Category")
            
            
            fig.suptitle(f'Proportion of Categories for {parameter_name} ({parameter_stable_id}), {sex}s', size = 26, ha = 'center',  y=1.001)
            plt.tight_layout()
            
            plt.savefig(f'{save}/barplot_genes_combined_{parameter_stable_id}_{sex}_rel_categories.png', dpi = 500, bbox_inches='tight')
            
            if save_memory == True:
                plt.clf() #clear current figure, save memory
        print("----------")
        print()
    print("-----------------------")
    print()
        
#%%
gene_symbols = ["Ap4e1","Dbn1","Nxn","Prkab1"]
plot_cat_combined(gene_symbols, sexes = ["female"], parameter_stable_ids = ["IMPC_EYE_017_001"])