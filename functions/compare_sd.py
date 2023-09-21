# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 12:27:59 2023

@author: Leonard Borchardt
"""

import os
import warnings
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import preprocessing as pre

#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')
#%%

def compare_sd(BW_age_in_weeks: int = 14, min_centers: int = 4, sd_calculated: bool = False, id_plots_done: bool = False, save_memory: bool = True):
    
    '''
    
    compares the measurement accuracy for each parameter between phenotyping centers
    by calculating and plotting (bar plot) the normalized center standard deviation for the recent control data of each center
    
    
    modules needed:
        import os
        import warnings
        import pandas as pd
        import scipy.stats as stats
        import matplotlib.pyplot as plt
        import seaborn as sns

    
    Parameters
    ----------
    BW_age_in_weeks: int, optional, relevant for IMPC_BWT_008_001, calculations for specific age [week], default: 14
    min_centers: int, min. phenotyping centers for each sd comparison, default: 4, optional
    sd_calculated: bool, skips calculation of group statistics if data frame already available, default: False, optional
    id_plots_done: bool, skips plotting of individual parameter stable id standard deviations, default: False, optional
    save_memory: bool, saves memory if "True" (plots get saved only), default: True, optional
    

    Returns (and saves):
    -------
    None
    automatically assigns a variable for a data frame with WT statistics for each phenotyping center
    
    additional saving of the data frame
    '''
    
    # ignore warnings
    warnings.filterwarnings("ignore")

    
    save = r'./overall_analysis'
    # check whether directory already exists, otherwise create new folder
    if not os.path.exists(save):
        os.mkdir(save)

    save = r'./overall_analysis/standard_deviation'
    # check whether directory already exists, otherwise create new folder
    if not os.path.exists(save):
        os.mkdir(save)



    if sd_calculated == False:

        read = r'./controls'
        
        
        ### compare standard deviations ###
        print("compare standard deviation between centers")
        print("------------------------------------------")
        print()
        
        controls = pd.read_csv(f'{read}/recent_controls.csv')
        
        #unidimensional parameters at most centers, + "IMPC_BWT_008_001", - "IMPC_BWT_001_001"
        ps_ids= [
        "IMPC_ABR_004_001","IMPC_ABR_006_001","IMPC_ABR_008_001","IMPC_ABR_010_001","IMPC_ABR_012_001","IMPC_ACS_001_001",
        "IMPC_ACS_002_001","IMPC_ACS_003_001","IMPC_ACS_004_001","IMPC_ACS_006_001","IMPC_ACS_007_001","IMPC_ACS_008_001",
        "IMPC_ACS_009_001","IMPC_ACS_033_001","IMPC_ACS_034_001","IMPC_ACS_035_001","IMPC_ACS_037_001","IMPC_CAL_001_001",
        "IMPC_CAL_002_001","IMPC_CAL_017_001","IMPC_CBC_004_001","IMPC_CBC_005_001","IMPC_CBC_006_001","IMPC_CBC_007_001",
        "IMPC_CBC_008_001","IMPC_CBC_009_001","IMPC_CBC_010_001","IMPC_CBC_012_001","IMPC_CBC_013_001","IMPC_CBC_014_001",
        "IMPC_CBC_015_001","IMPC_CBC_016_001","IMPC_CBC_017_001","IMPC_CBC_018_001","IMPC_CSD_032_001","IMPC_DXA_001_001",
        "IMPC_DXA_002_001","IMPC_DXA_003_001","IMPC_DXA_004_001","IMPC_DXA_005_001","IMPC_DXA_007_001","IMPC_DXA_008_001",
        "IMPC_DXA_009_001","IMPC_DXA_010_001","IMPC_GRS_003_001","IMPC_GRS_008_001","IMPC_GRS_009_001","IMPC_GRS_010_001",
        "IMPC_GRS_011_001","IMPC_HEM_001_001","IMPC_HEM_002_001","IMPC_HEM_003_001","IMPC_HEM_004_001","IMPC_HEM_005_001",
        "IMPC_HEM_006_001","IMPC_HEM_007_001","IMPC_HEM_008_001","IMPC_HWT_007_001","IMPC_HWT_008_001","IMPC_HWT_012_001",
        "IMPC_IPG_001_001","IMPC_IPG_010_001","IMPC_IPG_011_001","IMPC_IPG_012_001","IMPC_OFD_009_001","IMPC_OFD_010_001",
        "IMPC_OFD_012_001","IMPC_OFD_014_001","IMPC_OFD_016_001","IMPC_OFD_020_001","IMPC_PAT_049_002","IMPC_BWT_008_001"
        ]
        
        #only use the data for this ids for meaningful sd calculation (e.g. no series parameters incl.)
        controls = controls[controls["parameter_stable_id"].isin(ps_ids)]
        
        #drop nan values
        controls = controls.dropna(subset=['data_point'])
        
        control_stats_list = []
        
        for sex in ["male","female"]:
            
            controls_sex = controls[controls.sex == sex]
            
            
            # for IMPC_BWT_008_001, only consider mice of same age (default: week 14)
            controls_sex = controls_sex.drop(controls_sex[(controls_sex.parameter_stable_id == "IMPC_BWT_008_001") & (controls_sex.age_in_weeks != BW_age_in_weeks)].index)
            
      
            ##remove multiple measurements
            #by forming the mean of data_point for each specimen (---> 1 data point for each specimen)
            c = controls_sex.groupby(["sex","specimen_id","phenotyping_center","parameter_stable_id","parameter_name"])['data_point'].mean().to_frame().reset_index()
                #to_frame() bcs. output is a Series, reset_index() to not have data_point as index
                
                            
            c = c.sort_values(by=["parameter_stable_id",'phenotyping_center']).reset_index(drop=True)
            #print(controls)
            
        
            
            parameter_stable_ids = []
            parameter_names = []
            centers = []
            sample_sizes = []
            means = []
            medians = []
            sds = []
            Shapiro_tests = []
            CIs = []
            
            controls_gr = c.groupby(["parameter_stable_id","phenotyping_center"])
                
            for name, group in controls_gr:
                    
                # skip very small groups for meaningful statistics
                if len(group) < 3:
                    continue
                
                parameter_stable_ids.append(name[0])
                centers.append(name[1])
                parameter_names.append(list(group.parameter_name)[0])
            
            
                ###descriptive statistics
                
                #sample size
                sample_size = len(set(group.specimen_id))
                sample_sizes.append(sample_size)
                      
                #group mean
                mean = group.data_point.mean()
                means.append(mean)
                        
                #group median
                median = group.data_point.median()
                medians.append(median)
                        
                # calculate standard deviation for each group
                sd = group.data_point.std()
                sds.append(sd)
                
                #check for normal distribution
                W, p = stats.shapiro(group.data_point)
                Shapiro_tests.append(p)
                        
                # 95 % CI for the true mean of each group (CI for small sample sizes using t distribution)
                if len(group.data_point) < 30:
                    #CI sample sizes < 30 using t distribution
                    CI = stats.t.interval(confidence=0.95, df=len(group.data_point)-1, loc=group.data_point.mean(), scale=stats.sem(group.data_point))
                        #confidence instead of alpha!
                else:
                    CI = stats.norm.interval(confidence=0.95, loc=group.data_point.mean(), scale=stats.sem(group.data_point))
                CIs.append(CI)
    
                # birth time frame of the controls #
                controls_sex['date_of_birth'] = pd.to_datetime(controls_sex['date_of_birth'])
                       # Convert date_column to datetime format
                time_frame = controls_sex[(controls_sex.parameter_stable_id == name[0])&(controls_sex.phenotyping_center == name[1])].sort_values(by=["date_of_birth"]).reset_index(drop=True).date_of_birth[0], controls_sex[(controls_sex.parameter_stable_id == name[0])&(controls_sex.phenotyping_center == name[1])].sort_values(by=["date_of_birth"]).reset_index(drop=True).date_of_birth.iat[-1]
                    #tuple
        
                print(f'control data for {sex}, {name[0]} and {name[1]} successfully evaluated')
                
              
            control_stats = pd.DataFrame(dict(
                parameter_stable_id = parameter_stable_ids,
                parameter_name = parameter_names,
                sex = [sex]*len(parameter_stable_ids),
                phenotyping_center = centers,
                sample_size = sample_sizes,
                mean = means,
                median = medians,
                standard_deviation = sds,
                Shapiro_normality_p = Shapiro_tests,
                confidence_interval_true_mean = CIs,
                time_frame = [time_frame]*len(parameter_stable_ids)
                ))
            
            control_stats_list.append(control_stats)
            
        
        control_stats = pd.concat(control_stats_list, ignore_index=True)
        control_stats = control_stats.sort_values(by = ["parameter_stable_id", "sex","phenotyping_center"]).reset_index(drop=True)
        
        print('all control data successfully evaluated')  
        control_stats.to_csv(f"{save}/control_stats.csv", index=False)
    
    
    
    
    
    else:
        control_stats = pd.read_csv(f'{save}/control_stats.csv')
    
    
    
    ##########
    print("normalize (minmax) the standard deviation")
    
    #remove entries with sd = 0
    control_stats = control_stats.drop(control_stats[control_stats.standard_deviation==0].index)
    
    #normalize sd between 0 & 1 for each group
    stats_gr = control_stats.groupby(["parameter_stable_id","sex"])
    groups = []
    for name, group in stats_gr:
        # compare only groups with minimum number of centers
        if len(set(group.phenotyping_center)) < min_centers:
            continue
        sd_norm = pre.MinMaxScaler().fit_transform(group["standard_deviation"].array.reshape(-1, 1))
        group.insert(8, 'sd_normalized', sd_norm)
        groups.append(group)
    
    control_stats_norm = pd.concat(groups, ignore_index=True)
    
    control_stats_norm.to_csv(f"{save}/control_stats_normalized.csv", index=False)
    globals()["control_stats_normalized"] = control_stats_norm
    
    print("standard deviation normalized")
    
    
    
    if id_plots_done == False:
        #### bar plots for parameter_stable_id_specific normalized sd comparison between centers ####
        print("make bar charts for parameter_stable_id_specific normalized sd comparison between centers") 
        
        save = r'./overall_analysis/standard_deviation/plots_sd_norm'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(save):
            os.mkdir(save)
        
        for parameter_stable_id in set(control_stats_norm.parameter_stable_id):
            print(f"bar chart for {parameter_stable_id}")
            
            df = control_stats_norm[control_stats_norm.parameter_stable_id == parameter_stable_id]
            parameter_name = list(df.parameter_name)[0]
            
            plt.figure(figsize=(8, 8), dpi = 350)
        
            sns.set_theme()
        
            bar = sns.catplot(data=df, x="phenotyping_center", y="sd_normalized", col="sex", kind="bar")
        
            bar.fig.subplots_adjust(top=.85)
        
            if parameter_stable_id == "IMPC_BWT_008_001":
                bar.figure.suptitle(f'Normalized standard deviation across centers for {parameter_name} ({parameter_stable_id}), age: {BW_age_in_weeks} weeks', size = 18, ha = 'center')
                #bug: ha/horizontalaligment= 'right' ----> 'left'
        
            else:
                bar.figure.suptitle(f'Normalized standard deviation across centers for {parameter_name} ({parameter_stable_id})', size = 18, ha = 'center')
                #bug: ha/horizontalaligment= 'right' ----> 'left'
            
            for ax in bar.axes.flat:
                ax.tick_params(axis='x', labelrotation = 90)
                ax.set_xlabel("Phenotyping Center", fontsize=14)
                ax.set_ylabel("Normalized Standard Deviation", fontsize=14)
            
                #for i in bar.containers: ## similar to that
                    ###bar.bar_label(i,)
            
            if parameter_stable_id == "IMPC_BWT_008_001":
                
                plt.savefig(f"{save}/bar_plot_{parameter_stable_id}_normalized_standard_deviation_{BW_age_in_weeks}_weeks.png", dpi = 500, bbox_inches='tight')
                #whole path for plot saving important
                
            else:
                plt.savefig(f"{save}/bar_plot_{parameter_stable_id}_normalized_standard_deviation_.png", dpi = 500, bbox_inches='tight')
                #whole path for plot saving important
            
            if save_memory == True:
                plt.clf()
                #clear current figure, save memory
        
        
    print()    
        
    #### mean of normalized sd for each center and sex
    print("bar chart for sex-center-means of normalized standard deviations")
    
    save = r'./overall_analysis/standard_deviation'
    
    control_stats_normalized_sex_mean = control_stats_norm.groupby(["sex", "phenotyping_center"])["sd_normalized"].mean().to_frame().reset_index()
    control_stats_normalized_sex_mean.to_csv(f"{save}/control_stats_normalized_sex_mean.csv", index=False)
    globals()["control_stats_normalized_sex_mean"] = control_stats_normalized_sex_mean
    
    plt.figure(figsize=(8, 8), dpi = 350)

    sns.set_theme()

    bar = sns.catplot(data=control_stats_normalized_sex_mean, x="phenotyping_center", y="sd_normalized", col="sex", kind="bar")

    bar.fig.subplots_adjust(top=.85)

    bar.figure.suptitle('Normalized standard deviation sex mean over all parameters', size = 18, ha = 'center')
        #bug: ha/horizontalaligment= 'right' ----> 'left'
    
    for ax in bar.axes.flat:
        ax.tick_params(axis='x', labelrotation = 90)
        ax.set_xlabel("Phenotyping Center", fontsize=14)
        ax.set_ylabel("Normalized Standard Deviation", fontsize=14)
    
        #for i in bar.containers: ## similar to that
            ###bar.bar_label(i,)
    
    plt.savefig(f"{save}/bar_plot_normalized_standard_deviation_sex_mean.png", dpi = 500, bbox_inches='tight')
    #whole path for plot saving important
     
        
     
    #### mean of normalized sd for each center
    print("bar chart for center-means of normalized standard deviations")
    
    save = r'./overall_analysis/standard_deviation'
    
    control_stats_norm_mean = control_stats_norm.groupby(["phenotyping_center"])["sd_normalized"].mean().to_frame().reset_index()
    control_stats_norm_mean.to_csv(f"{save}/control_stats_normalized_mean.csv", index=False)
    globals()["control_stats_normalized_mean"] = control_stats_norm_mean
    
    
    plt.figure(figsize=(8, 14), dpi = 350)

    sns.set_theme()

    bar = sns.catplot(data=control_stats_norm_mean, x="phenotyping_center", y="sd_normalized", kind="bar")

    bar.fig.subplots_adjust(top=.85)

    bar.figure.suptitle('Normalized standard deviation mean over all parameters', size = 18, ha = 'center')
        #bug: ha/horizontalaligment= 'right' ----> 'left'
    
    for ax in bar.axes.flat:
        ax.tick_params(axis='x', labelrotation = 90)
        ax.set_xlabel("Phenotyping Center", fontsize=14)
        ax.set_ylabel("Normalized Standard Deviation", fontsize=14)
    
        #for i in bar.containers: ## similar to that
            ###bar.bar_label(i,)
    
    plt.savefig(f"{save}/bar_plot_normalized_standard_deviation_mean.png", dpi = 500, bbox_inches='tight')
    
#%%
compare_sd()