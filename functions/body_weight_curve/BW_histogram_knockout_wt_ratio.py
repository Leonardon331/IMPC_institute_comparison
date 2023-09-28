# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 12:19:41 2023

@author: Leonard Borchardt
"""

import os
import warnings
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import math

#set wd:
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')

#%%

def BW_histogram_knockout_wt_ratio(gene_symbols: list, parameter_stable_id: str = "IMPC_BWT_008_001", focus_zygosity: bool = False, less_outliers: bool = False, min_week: int = "None", max_week: int = "None", save_memory = True):
    '''
    needs single parameter gene dataframe (.csv) with knockout_wt_ratio_column, has to be in a parameters/ps_id/gene folder in the cwd (current working directory)

    description:
    reads in single parameter gene dataframe (.csv) with knockout_wt_ratio_column in cwd/parameters/id/gene and plots/saves kwr histograms for sex and zygosity across center
    
    BW_...: aditional feature: filters after age greater/equal min_week and smaller/equal max _week parameter to make the results less fluctuating (in early age)
    
    modules needed:
     - import os
     - import warnings
     - import pandas as pd
     - import scipy.stats as stats
     - import matplotlib.pyplot as plt
     - from matplotlib.ticker import MaxNLocator
     - import seaborn as sns
     - import math

    Parameters
    ----------
    gene_symbols: list of strings, symbols of the genes used for center comparison
    parameter_stable_id: str, IMPC pipline parameter id that is connected to body weight curve obtaining, default: "IMPC_BWT_008_001", id allows parameter measurement comparison between centers
    focus_zygosity: bool, default: False; if "True", comparison focuses on zygosity with most data across centers
    less_outliers: remove outliers 3xIQR from border quartiles, False by default, optional
    min_week: int, min. mice age, default: "None", optional
    max_week: int, max. mice age, default: "None", optional
    save_memory: saves memory if "True" (plots get saved only), default: True, optional
    
    Returns
    -------
    None

    '''
    
    # ignore warnings
    warnings.filterwarnings("ignore")
    
    for gene_symbol in gene_symbols:
        
        print()
        print(gene_symbol)
        print(parameter_stable_id)
        
        
        
        ########
        #create folder if it doesn't already exist
        save = r'./parameters'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(save):
            os.mkdir(save)
        
        #create subfolder if it doesn't already exist
        save = f'{save}/{parameter_stable_id}'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(save):
            os.mkdir(save)
        
        #create subfolder if it doesn't already exist
        save = f'{save}/{gene_symbol}'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(save):
            os.mkdir(save)
            
        #create subfolder if it doesn't already exist
        save = f'{save}/histograms'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(save):
            os.mkdir(save)
            
        #########
        
        
        
        #path for reading
        path = f'./parameters/{parameter_stable_id}/{gene_symbol}'
        
        
        #read in the data in the cwd
        try:
            df = pd.read_csv(f"{path}/{gene_symbol}_{parameter_stable_id}_kwr.csv")
        except FileNotFoundError:
            print("missing data for center comparison")
            continue
        
        #drop nan values
        df = df.dropna(subset=['knockout_wt_ratio'])
        
        parameter_name = df.parameter_name[0]
        
        ##remove multiple measurements
        #by forming the mean of data_point for each specimen age in weeks (---> 1 data point for each specimen age in weeks)
        df = df.groupby(["specimen_id","phenotyping_center","sex","zygosity","age_in_weeks"])['knockout_wt_ratio'].mean().to_frame().reset_index()
            #to_frame() bcs. output is a Series, reset_index() to not have knockout_wt_ratio as index    
     
        
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
        
        
        ##### optional filtering: min_week
        #filter after min age in weeks (if min_week parameter in use)
        if min_week != "none":
            df = df[df.age_in_weeks >= min_week]
        #####
        
        ##### optional filtering: max_week
        #filter after max age in weeks (if max_week parameter in use)
        if max_week != "none":
            df = df[df.age_in_weeks <= max_week]
        #####
        
        
        
        ##less outliers, default = False
        #by removing points with > 3x IQR distance from lower/upper quartile
        if less_outliers == True:
            groups = []
            df_groups = df.groupby(["phenotyping_center", "sex", "zygosity"])
            for name,group in df_groups:
                #(without name: interation variable equals tuples) 
                #print(group.knockout_wt_ratio)
                #calculate interquartile range (IQR) and drop measurements that are
                #more than 3x higher/lower the first/third quartile
                iqr = stats.iqr(group.knockout_wt_ratio)
                first_quartile = group.knockout_wt_ratio.quantile(.25)
                third_quartile = group.knockout_wt_ratio.quantile(.75)
                outer_fence = [first_quartile - 3*iqr, third_quartile + 3*iqr]
                for specimen in group.specimen_id:
                    for week in group.age_in_weeks:
                        #print(group[(group.specimen_id == specimen)&(group.age_in_weeks == week)].knockout_wt_ratio)
                        try:
                            if group[(group.specimen_id == specimen)&(group.age_in_weeks == week)].knockout_wt_ratio.item() < outer_fence[0] or group[(group.specimen_id == specimen)&(group.age_in_weeks == week)].knockout_wt_ratio.item() > outer_fence[1]:
                                #check whenever this data point from this specimen/week lays outside the boundries
                                group.drop(group[(group.specimen_id==specimen)&(group.age_in_weeks==week)].index, inplace = True)
                        except ValueError:
                            continue
                            
                groups.append(group)
                
            df = pd.concat(groups,ignore_index=True)
           
                
        #remove groups with less than 3 entries for meaningful statistics
        groups = []
        df_groups = df.groupby(["phenotyping_center", "sex", "zygosity"])
        for name,group in df_groups:
            if len(group) < 3:
                del group
            else:
                groups.append(group)   
                
        df = pd.concat(groups,ignore_index=True)
            
        df = df.sort_values(by=['zygosity', 'sex', 'phenotyping_center']).reset_index(drop=True)
        #print(df)
        
        
        ##### histogram ####
        
        for sex in set(df.sex):
            for zygosity in set(df[df.sex == sex]["zygosity"]): 
                
                df_sz = df[(df.sex == sex)&(df.zygosity == zygosity)]
                center_count = len(set(df_sz["phenotyping_center"]))
                
                centers = list(set(df_sz.phenotyping_center))
                centers.sort()
                
                # Increase the height of each subplot
                subplot_height = 4  # You can adjust this value to control the height
                subplot_count_per_row = 3
                subplot_rows = math.ceil(center_count / subplot_count_per_row)
                total_fig_height = subplot_height * subplot_rows
                
                # Create subplots with adjusted figsize
                fig, axs = plt.subplots(subplot_rows, subplot_count_per_row, figsize=(12, total_fig_height), dpi=500)
                
                sns.set_theme()
            
                if center_count == 1:
                    fig.suptitle(f'{gene_symbol}/{parameter_name} ({parameter_stable_id})/{zygosity}/{sex}; week: min = {min_week}, max = {max_week}', size = 15, ha = 'right')
                    #bug: ha/horizontalaligment= 'right' ----> 'left'yout()
                else:
                    fig.suptitle(f'{gene_symbol}/{parameter_name} ({parameter_stable_id})/{zygosity}/{sex}; week: min = {min_week}, max = {max_week}', size = 15, ha = 'center')
                
                i = 0
                for ax in axs.flat:
                    #.flat creates 1 D array
                    if i >= center_count:
                        ax.remove()
                        i += 1
                        continue
                
                    d = df_sz[df_sz.phenotyping_center == centers[i]] 
                        #data for certain sex/zygosity/phenotyping center
                    sns.histplot(data=d, x="knockout_wt_ratio", kde=True, ax=ax)
                    
                    ax.set_title(f'{centers[i]}')
                    
                    #ax.tick_params(axis='x', labelrotation = 45)
                    ax.set_xlabel("Knockout-WT-Ratio", fontsize=10)
                    ax.set_ylabel("Count", fontsize=10)
                    #ax.label_outer()
                        #labels only for outer axis
                    
                    # Set y-axis ticks as integers
                    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
                    
                    i += 1
                    
                plt.tight_layout()
                
                if less_outliers == True:
                    plt.savefig(f"{save}/histogram_{gene_symbol}_{parameter_stable_id}_kwr_{zygosity}_{sex}_week_min_{min_week}_max_{max_week}.png", dpi = 500, bbox_inches='tight')
                    #whole path for plot saving important
                else:
                    plt.savefig(f"{save}/histogram_{gene_symbol}_{parameter_stable_id}_kwr_{zygosity}_{sex}_week_min_{min_week}_max_{max_week}_outliers_incl.png", dpi = 500, bbox_inches='tight')
                
                
                
                if save_memory == True:
                    plt.clf()
                    #clear current figure
#%%

genes = ["Prkab1", "Dbn1", "Ap4e1", "Nxn"]
BW_histogram_knockout_wt_ratio(genes, min_week = 14, max_week = 16)