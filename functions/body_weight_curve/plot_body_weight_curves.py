# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 12:02:23 2023

@author: Leonard Borchardt
"""
import os
import warnings
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')
#%%

def plot_body_weight_curves(gene_symbols: list, parameter_stable_id: str = "IMPC_BWT_008_001", focus_zygosity: bool = True, mean_only: bool = True):
    '''
    reads in gene data frame in genes folder in cwd and
    plots body weight curves (absolute values) for given genes for each sex for most measured zygosity across centers rounded to the nearest week and
    saves the plot
    
    import modules:
      import os
      import warnings
      import pandas as pd
      import seaborn as sns  
      import matplotlib.pyplot as plt 
    
    parameters:
     - gene_symbol: list of strings, symbols for genes
     - parameter_stable_id: str, default: "IMPC_BWT_008_001" (IMPC parameter ID for body weight curve obtaining)
     - focus_zygosity: bool, default: True, focuses comparison on zygosity with most data across centers, optional
     - mean_only: bool, plot only 1 point per week for mean of each week/zygosity/sex/phenotyping center, True by default, optional
     
    returns:
     - None
    '''
    
    # ignore warnings
    warnings.filterwarnings("ignore")
    
    print(parameter_stable_id)
    
    for gene_symbol in gene_symbols:
        print(gene_symbol)
        
        ########
        #create folder if it doesn't already exist
        path = r'./overall_analysis'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(path):
            os.mkdir(path)
        
        #create subfolder for saving if it doesn't already exist
        path = f'{path}/body_weight_curves'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(path):
            os.mkdir(path)
        #########
        
        #read in and filter gene data    
        df =  pd.read_csv(f"./genes/{gene_symbol}/{gene_symbol}.csv")
        df = df[(df.parameter_stable_id == parameter_stable_id)]
        
        # make measured values positive
        df['data_point'] = df['data_point'].abs()
        
        
        #drop nan values
        df = df.dropna(subset=['data_point'])
        
        #remove multiple measurements
        #by forming the mean of data_point for each specimen age in weeks (---> 1 data point for each specimen age in weeks)
        df = df.groupby(["specimen_id","phenotyping_center","sex","zygosity", "age_in_weeks"])['data_point'].mean().to_frame().reset_index()
            #to_frame() bcs. output is a Series, reset_index() to not have data_point as index    
        
        
        
        #### additional filtering if parameter focus_zygosity = True (default)
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
          
        if focus_zygosity == False:
            zygosity = "all_zygosities"
            

        if mean_only == True:
        #calculate the mean for each week/center/zygosity/sex
            
            ##first: less strong outliers for better mean
            #by removing points with > 3x IQR distance from lower/upper quartile
            groups = []
            df_groups = df.groupby(["phenotyping_center", "sex", "zygosity", "age_in_weeks"])
            for name,group in df_groups:
                #(without name: interation variable equals tuples) 
                #print(group.data_point)
                #calculate interquartile range (IQR) and drop measurements that are
                #more than 3x higher/lower the first/third quartile
                iqr = stats.iqr(group.data_point)
                first_quartile = group.data_point.quantile(.25)
                third_quartile = group.data_point.quantile(.75)
                outer_fence = [first_quartile - 3*iqr, third_quartile + 3*iqr]
                for specimen in group.specimen_id:
                    if group[group.specimen_id == specimen].data_point.item() < outer_fence[0] or group[group.specimen_id == specimen].data_point.item() > outer_fence[1]:
                        group.drop(group[group.specimen_id == specimen].index, inplace = True)
            
                groups.append(group)
            
            df = pd.concat(groups,ignore_index=True)
            
            
            #now calculate mean
            df = df.groupby(["phenotyping_center", "sex", "zygosity", "age_in_weeks"])["data_point"].mean().to_frame().reset_index()
                #now df only with week means for each center/zygosity/sex
    
    
    
        ### plotting
        
        plt.figure(figsize=(14, 14), dpi = 500)
    
        sns.set_theme()
        
        #age rounded to the neares week
        rel = sns.relplot(
            data=df, x="age_in_weeks", y="data_point", col="sex", row = "zygosity", 
            hue="phenotyping_center", kind="line", errorbar = "sd")
        
        #rename the legend title
        rel._legend.set_title("Phenotyping Center")
        
        #move the following title higher
        if len(set(df.zygosity)) == 1:
            rel.fig.subplots_adjust(top=.85)
        else:
            rel.fig.subplots_adjust(top=.92)
    
        rel.fig.suptitle(f'Body weight curve ({parameter_stable_id}): {gene_symbol}', size = 18, ha = 'right')
        #bug: ha/horizontalaligment= 'right' ----> 'left'
        
        for ax in rel.axes.flat:
            ax.set_xlabel("Age [W]", fontsize=10)
            ax.set_ylabel("Body weight [g]", fontsize=10)
        
        
        
        if mean_only == True:
            plt.savefig(f"{path}/body_weight_curves_{gene_symbol}_{parameter_stable_id}_{zygosity}.png", dpi = 500, bbox_inches='tight')
            
        if mean_only == False:
            plt.savefig(f"{path}/body_weight_curves_{gene_symbol}_{parameter_stable_id}_{zygosity}_sd.png", dpi = 500, bbox_inches='tight')
#%%
plot_body_weight_curves(["Prkab1", "Dbn1", "Ap4e1", "Nxn"], focus_zygosity=False)

