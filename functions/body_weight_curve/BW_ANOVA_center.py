# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 10:17:32 2023

@author: Leonard Borchardt
"""

#one-way ANOVA for knockout_wt_ratio with factor phenotyping center


import os
import warnings
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from tabulate import tabulate
import pingouin as pg
import scikit_posthocs as sp


#set wd:
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')
#%%
def BW_ANOVA_center(gene_symbols: list, parameter_stable_id: str = "IMPC_BWT_008_001", post_hoc_method: str = "ttest", p_adjust: str = "fdr_bh", plot_type: str = "box", focus_zygosity: bool = True, less_outliers: bool = False, min_week: int = "None", max_week: int = "None", save_memory: bool = False):
    '''
    needs single parameter gene dataframe (.csv) with knockout_wt_ratio_column, has to be in a parameters/ps_id/gene folder in cwd (current working directory)

    description:
    reads in single parameter gene dataframe (.csv) with knockout_wt_ratio_column in cwd and prints a table & returns a data frame for F-statistic and p-value
    for one way ANOVA with independent variable/factor "phenotyping_center" and dependent variable knockout_wt_ratio, the returned dataframe contains additional
    group-descriptive statistic parameters
    
    BW_...: aditional feature: filters after age greater/equal min_week and smaller/equal max _week parameter to make the results less fluctuating (in early age)
    
    modules needed:
     - import os
     - import warnings
     - import pandas as pd
     - import scipy.stats as stats
     - import matplotlib.pyplot as plt
     - import seaborn as sns
     - from tabulate import tabulate
     - import pingouin as pg
     - import scikit_posthocs as sp

    Parameters
    ----------
    gene_symbols: list of strings, symbols of the genes used for center comparison
    parameter_stable_id: str, IMPC pipline parameter id that is connected to body weight curve obtaining, default: "IMPC_BWT_008_001", id allows parameter measurement comparison between centers, optional
    post_hoc_method: str, "ttest" (default), perform multi Welch’s t-test
    
    p_adjust: str, method for adjusting the post hoc p value (dealing with inflated Type I error for multi tests),
                Available methods are: None, ‘bonferroni’: one-step correction ‘sidak’: one-step correction ‘holm-sidak’: step-down method using Sidak adjustments 
                ‘holm’: step-down method using Bonferroni adjustments ‘simes-hochberg’ : step-up method (independent) ‘hommel’: closed method based on Simes tests (non-negative)
                ‘fdr_bh’: Benjamini/Hochberg (non-negative) ‘fdr_by’: Benjamini/Yekutieli (negative) ‘fdr_tsbh’: two stage fdr correction (non-negative)
                ‘fdr_tsbky’: two stage fdr correction (non-negative), optional
                default: "fdr_bh"
    
    plot_type: str, categorical plot type, default: "box"; other options: "violin", "swarm", "strip", "boxen", optional
    focus_zygosity: bool, comparison focuses on zygosity with most data across centers, default: True, optional
    less_outliers: bool, remove outliers 3xIQR from border quartiles, default: False, optional
    min_week: int, min. mice age, default: "None", optional
    max_week: int, max. mice age, default: "None", optional
    save_memory: bool, saves memory (no returns only saving), default: False, optional
    
    Returns
    -------
    None
    assignes automatically variable for dataframe with ANOVA_results and group-descriptive statistics and for significant ANOVAs also post hoc data frames 
    
    prints:
    table for F-statistic and p-value
    for one way ANOVA with independent variable/factor "phenotyping_center" and dependent variable knockout_wt_ratio

    '''
    
    # ignore warnings
    warnings.filterwarnings("ignore")
    
    for gene_symbol in gene_symbols:
        
        print()
        print()
        print()
        print(gene_symbol)
        print(parameter_stable_id)
        print()
        
        #path for reading and saving
        path = f'./parameters/{parameter_stable_id}/{gene_symbol}'
        
        
        #read in the data in the cwd
        try:
            df = pd.read_csv(f"{path}/{gene_symbol}_{parameter_stable_id}_kwr.csv")
        except FileNotFoundError:
            print("missing knockout-wt-ratio data, run function compare_knockout_to_wt.py first")
            continue
        
        parameter_name = df.parameter_name[0]
        
        #check if "knockout_wt_ratio" column contains data
        if df['knockout_wt_ratio'].empty or df['knockout_wt_ratio'].isna().all():
            print(f"{parameter_stable_id} seems to be non-numerical")
            print(f"skipping {parameter_stable_id}")
            continue
        
        #drop nan values
        df = df.dropna(subset=['knockout_wt_ratio'])
        
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
        
        
        
        ##less outliers, default: False
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
           
                
        #remove groups with less than 3 samples for meaningful statistics
        groups = []
        df_groups = df.groupby(["phenotyping_center","sex", "zygosity"])
        for name,group in df_groups:
            if len(group) < 3: #sample size check
                continue
            else:
                groups.append(group) 
                
        df = pd.concat(groups,ignore_index=True)
        
        #remove groups with less than 2 centers for meaningful comparison
        groups = []
        df_groups = df.groupby(["sex", "zygosity"])
        for name,group in df_groups:
            if len(set(group.phenotyping_center)) < 2: #center number check
                continue
            else:
                groups.append(group) 
                
        df = pd.concat(groups,ignore_index=True)
            
            
          
        df = df.sort_values(by=['zygosity', 'sex', 'phenotyping_center']).reset_index(drop=True)
        #print(df) 
        
        
        
        ###### iterate through every sex and zygosity and compare them across centers
        results = []
        all_post_hoc_p = []
        
        for sex in set(df.sex):
            for zygosity in set(df[df.sex == sex]["zygosity"]):
                
                d = df[(df.sex == sex) & (df.zygosity == zygosity)]
                
                print()
                print(f"{zygosity} {sex} knockout_wt_ratio ANOVA across phenotyping centers")
                print()
                print("check test criterias first and raise warning if not met")
                
                Shapiro_tests = []
                print()
                print(f"Shapiro normality tests for each {zygosity} {sex} phenotyping center sample group")
                
                for center in set(d["phenotyping_center"]):
                    
                    try:
                        W,p = stats.shapiro(d[d.phenotyping_center == center].knockout_wt_ratio)
                        #if p > 0.05:
                            #print(f"knockout_wt_ratio for {zygosity} {sex}s from {center} probably normaly distributed")
                            #print()  
                        if p < 0.05:
                            print()
                            print(f"knockout_wt_ratio for {zygosity} {sex}s from {center} may not be normally distributed")
                            print(f"W-statistic: {W}")
                            print(f"p-value: {p}")
                        
                    except ValueError:
                        print("sample size too small")
                        p = "sample size too small"
                    
                    Shapiro_tests.append(p)
                
                #print(f"check {zygosity} {sex}s across phenotyping centers for variance homogeneity:")
                #check for variance homogeneity
                centers_sex_zygosity = d.groupby('phenotyping_center')['knockout_wt_ratio']
                #groups of knockout_wt_ratio Series
                
                ## Perform Levene's test for variance homogeneity
                #first get individual center groups
                knockout_wt_ratios = []
                #list of arrays
                for gr_name, gr_df in centers_sex_zygosity:
                    # extract kwr values of group (center) into array
                    #(by getting df for groups and the kwr values then as an array)
                    ##print(group_df)
                    ##return group_df
                    knockout_wt_ratios.append(gr_df.to_numpy())
                
                try:
                    levene_statistic, levene_p_value = stats.levene(*knockout_wt_ratios)
                    #unpack list of arrays
                    #print(f"Levene's test statistic: {levene_statistic}")
                    #print(f"p-value: {levene_p_value}")
                    #if levene_p_value > 0.05:
                        #print(f"variances between centers for {zygosity} {sex}s probably equal")
                        #print()
                        
                    #levene warning (not necessary)
                    #if levene_p_value < 0.05:
                        #print(f"variances between centers for {zygosity} {sex}s may not be equal")
                        #print(f"Levene's test statistic: {levene_statistic}")
                        #print(f"p-value: {levene_p_value}")
                        #print()
                
                except ValueError:
                    #print("not enough groups for ANOVA (need at least 2)")
                    levene_p_value = "not enough groups"
                
                
                
                
                ###descriptive statistics
                
                #sample size
                sample_sizes = []
                for center in set(d["phenotyping_center"]):
                    sample_size = len(centers_sex_zygosity.get_group(center))
                    sample_sizes.append(sample_size)
                
                
                #check group for normal distribution
                #Shapiro_tests = []
                #for center in set(d["phenotyping_center"]):
                    #center_values = centers_sex_zygosity.get_group(center)
                    #W, p = stats.shapiro(center_values)
                    #Shapiro_tests.append(p)
                
                
                #group mean
                means = []
                for center in set(d["phenotyping_center"]):
                    mean = centers_sex_zygosity.get_group(center).mean()
                    means.append(mean)
                    
                #group median
                medians = []
                for center in set(d["phenotyping_center"]):
                    median = centers_sex_zygosity.get_group(center).median()
                    medians.append(median)
                    
                # calculate standard deviation for each group
                sds = []
                for center in set(d["phenotyping_center"]):
                    sd = centers_sex_zygosity.get_group(center).std()
                    sds.append(sd)
                    
                # 95 % CI for the true mean of each group (CI for small sample sizes using t distribution)
                CIs = []
                for center in set(d["phenotyping_center"]):
                    center_values = centers_sex_zygosity.get_group(center)
                    #print(center_values)
                    #print(type(center_values))
                    if len(center_values) < 30:
                        #CI sample sizes < 30 using t distribution
                        CI = stats.t.interval(confidence=0.95, df=len(center_values)-1, loc=center_values.mean(), scale=stats.sem(center_values))
                            #confidence instead of alpha!
                    else:
                        CI = stats.norm.interval(confidence=0.95, loc=center_values.mean(), scale=stats.sem(center_values))
                        
                    CIs.append(CI)
                
                
                
                
                ###ANOVA
                print()
                print(f"perform ANOVA across centers for {zygosity} {sex} knockout_wt_ratio")
                
                try:
                    if levene_p_value > 0.05:
                    # homogeneous variance
                        print("homogene variance --> classic ANOVA")
                        ANOVA = pg.anova(data=d,dv="knockout_wt_ratio",between="phenotyping_center")
                        #print(ANOVA)
                        f_statistic = ANOVA.F.item()
                        #print(f_statistic)
                        p_value = ANOVA["p-unc"].item()
                        #print(p_value)
                        partial_eta_square = ANOVA["np2"].item()
                        
                    
                    if levene_p_value < 0.05:
                    # inhomogeneous variance
                        print("inhomogene variance --> Welch's ANOVA")
                        ANOVA = pg.welch_anova(data=d,dv="knockout_wt_ratio",between="phenotyping_center")
                        #print(ANOVA)
                        f_statistic = ANOVA.F.item()
                        #print(f_statistic)
                        p_value = ANOVA["p-unc"].item()
                        #print(p_value)
                        partial_eta_square = ANOVA["np2"].item()
                        
                    
                    result = pd.DataFrame({ 'gene_symbol': gene_symbol,
                                            'parameter_stable_id': parameter_stable_id,
                                            'parameter_name': parameter_name,
                                            'sex': sex,
                                            'zygosity': zygosity,
                                            'phenotyping_center': list(set(d.phenotyping_center)),
                                            'sample_size': sample_sizes,
                                            "Shapiro_normality_p": Shapiro_tests,
                                            "Levene_variance_homogeneity_p": levene_p_value,
                                            'mean': means,
                                            'median': medians,
                                            'standard_deviation': sds,
                                            'confidence_interval_true_mean': CIs,
                                            'F_Statistic_ANOVA': f_statistic,
                                            'p_Value_ANOVA': p_value,
                                            'partial_eta_square_ANOVA': partial_eta_square
                                            })
                

    
                except (ValueError, TypeError):
                    print(f'not enough data for ANOVA for {zygosity} {sex}s')
                        
                    result = pd.DataFrame({ 'gene_symbol': gene_symbol,
                                            'parameter_stable_id': parameter_stable_id,
                                            'parameter_name': parameter_name,
                                            'sex': sex,
                                            'zygosity': zygosity,
                                            'phenotyping_center': list(set(d.phenotyping_center)),
                                            'sample_size': sample_sizes,
                                            "Shapiro_normality_p": Shapiro_tests,
                                            "Levene_variance_homogeneity_p": levene_p_value,
                                            'mean': means,
                                            'median': medians,
                                            'standard_deviation': sds,
                                            'confidence_interval_true_mean': CIs,
                                            'F_Statistic_ANOVA': "not enough data",
                                            'p_Value_ANOVA': "not enough data",
                                            'partial_eta_square_ANOVA': "not enough data"
                                            })
                
                results.append(result)
                
                
                #### post hoc tests ###
                if p_value < 0.05:
                    print()
                    print("ANOVA significant")
                    
                    try:
                        print("perform post hoc tests between centers")
     
                        if post_hoc_method == "ttest":
                            print(f"method: multiple Welch’s t-tests with {p_adjust} p correction")
                            #using "holm" method to adjust the p value for multiple comparisons (keep Type I error low)
                            post_hoc_p = sp.posthoc_ttest(a = d, val_col = "knockout_wt_ratio", group_col = "phenotyping_center", equal_var = False, p_adjust = p_adjust, sort = True)
                                #perform Welch’s t-test, which does not assume equal population variance
                        
                        phenotyping_centers = list(set(d.phenotyping_center))
                        phenotyping_centers.sort()
                        
                        post_hoc_p.insert(0, "phenotyping_center", phenotyping_centers)
                        post_hoc_p.insert(0, "zygosity", zygosity)
                        post_hoc_p.insert(0, "sex", sex)
                        post_hoc_p.insert(0, "parameter_name", parameter_name)
                        post_hoc_p.insert(0, "parameter_stable_id", parameter_stable_id)
                        post_hoc_p.insert(0, "gene_symbol", gene_symbol)
                        #don't assign the data frame variable again (leads to error)
                        
                        all_post_hoc_p.append(post_hoc_p)
                        
                    except ZeroDivisionError:
                        continue
                print()
        
        try:
            results = pd.concat(results, ignore_index=True)
            #concat list out of data frames
            results = results.sort_values(by=['zygosity','sex','phenotyping_center']).reset_index(drop=True)
        except ValueError:
            print(f"no data to compare across centers for {gene_symbol}, {parameter_stable_id}, {sex}, {zygosity}")
            print(f'run filter_genes function with "{parameter_stable_id}" and check data_point column')
            continue
            #next gene
         
        
        if len(all_post_hoc_p) != 0:
            all_post_hoc_p = pd.concat(all_post_hoc_p, ignore_index=True)
    
        
        
        print()
        print("Interpretation of ANOVA results")
        
        #### create variables to access the ANOVA results
        
        homo_m = results[(results.sex == 'male') & (results.zygosity == 'homozygote')]
        if not homo_m.empty:
            try:
                homo_m_f = results[(results.sex == 'male') & (results.zygosity == 'homozygote')].F_Statistic_ANOVA.reset_index(drop=True)[0].item()
                    #.item() to convert Series object to scalar for tabulate funtion
                homo_m_p  = results[(results.sex == 'male') & (results.zygosity == 'homozygote')].p_Value_ANOVA.reset_index(drop=True)[0].item()
            except (TypeError, AttributeError):
                #cannot convert to scalar bcs of string character of entry or if [0] is already putting out a scalar (...why?)
                try:
                    homo_m_f = results[(results.sex == 'male') & (results.zygosity == 'homozygote')].F_Statistic_ANOVA.reset_index(drop=True)[0]
                    homo_m_p  = results[(results.sex == 'male') & (results.zygosity == 'homozygote')].p_Value_ANOVA.reset_index(drop=True)[0]
                except KeyError:
                    #no entry at all (bcs no data for particular sex and zygosity at all)
                    homo_m_f = "no data"
                    homo_m_p = "no data"
        else:
            homo_m_f = "no data"
            homo_m_p = "no data"
        
            
        homo_f = results[(results.sex == 'female') & (results.zygosity == 'homozygote')]
        if not homo_f.empty:
            try:
                homo_f_f = results[(results.sex == 'female') & (results.zygosity == 'homozygote')].F_Statistic_ANOVA.reset_index(drop=True)[0].item()
                homo_f_p = results[(results.sex == 'female') & (results.zygosity == 'homozygote')].p_Value_ANOVA.reset_index(drop=True)[0].item()
            except (TypeError, AttributeError):
                try:
                    homo_f_f = results[(results.sex == 'male') & (results.zygosity == 'homozygote')].F_Statistic_ANOVA.reset_index(drop=True)[0]
                    homo_f_p = results[(results.sex == 'female') & (results.zygosity == 'homozygote')].p_Value_ANOVA.reset_index(drop=True)[0]
                except KeyError:
                    homo_f_f = "no data"
                    homo_f_p = "no data"
        else:
            homo_f_f = "no data"
            homo_f_p = "no data"
        
            
        hetero_m = results[(results.sex == 'male') & (results.zygosity == 'heterozygote')]
        if not hetero_m.empty:
            try:
                hetero_m_f = results[(results.sex == 'male') & (results.zygosity == 'heterozygote')].F_Statistic_ANOVA.reset_index(drop=True)[0].item()
                hetero_m_p = results[(results.sex == 'male') & (results.zygosity == 'heterozygote')].p_Value_ANOVA.reset_index(drop=True)[0].item()
            except (TypeError, AttributeError):
                try:
                    hetero_m_f = results[(results.sex == 'male') & (results.zygosity == 'homozygote')].F_Statistic_ANOVA.reset_index(drop=True)[0]     
                    hetero_m_p = results[(results.sex == 'male') & (results.zygosity == 'heterozygote')].p_Value_ANOVA.reset_index(drop=True)[0]
                except KeyError:
                    hetero_m_f = "no data"
                    hetero_m_p = "no data"
        else:
            hetero_m_f = "no data"
            hetero_m_p = "no data"
            
        
        hetero_f = results[(results.sex == 'female') & (results.zygosity == 'heterozygote')]
        if not hetero_f.empty:
            try:
                hetero_f_f = results[(results.sex == 'female') & (results.zygosity == 'heterozygote')].F_Statistic_ANOVA.reset_index(drop=True)[0].item()
                hetero_f_p = results[(results.sex == 'female') & (results.zygosity == 'heterozygote')].p_Value_ANOVA.reset_index(drop=True)[0].item()
            except (TypeError, AttributeError):
                try:
                    hetero_f_f = results[(results.sex == 'male') & (results.zygosity == 'homozygote')].F_Statistic_ANOVA.reset_index(drop=True)[0]   
                    hetero_f_p = results[(results.sex == 'female') & (results.zygosity == 'heterozygote')].p_Value_ANOVA.reset_index(drop=True)[0]
                except KeyError:
                    hetero_f_f = "no data"
                    hetero_f_p = "no data"
        else:
            hetero_f_f = "no data"
            hetero_f_p = "no data"
        
            
        ### print the results in table format ###
        data = [["homozygote males between centers", homo_m_f, homo_m_p],
                ["homozygote females between centers", homo_f_f, homo_f_p], 
                ["heterozygote males between centers", hetero_m_f, hetero_m_p],
                ["heterozygote females between centers", hetero_f_f, hetero_f_p]]
      
        #define header names
        col_names = ["Groups", "F Statistic", "p Value"]
        #display table
        print(tabulate(data, headers=col_names, tablefmt="fancy_grid"))
        ###
        
        
        ####Interpretation of ANOVA results
        
        #ignoring missing data beween centers (not comparable)
        if isinstance(homo_m_p, str):
                homo_m_p = 1
        
        if isinstance(homo_f_p, str):
                homo_f_p = 1
        
        if isinstance(hetero_m_p, str):
                hetero_m_p = 1
        
        if isinstance(hetero_f_p, str):
                hetero_f_p = 1
    
    
        #making statements to center similarity
        if homo_m_p > 0.05 and homo_f_p > 0.05 and hetero_m_p > 0.05 and hetero_f_p > 0.05:
            print(f"the data generated for {gene_symbol} for parameter {parameter_stable_id} is across centers highly similar/comparable")
        
        if homo_m_p < 0.05:
            print(f"the data generated for {gene_symbol} for parameter {parameter_stable_id} is between centers for homozygote males not very similar/comparable")
            
        if homo_f_p < 0.05:
           print(f"the data generated for {gene_symbol} for parameter {parameter_stable_id} is between centers for homozygote females not very similar/comparable") 
        
        if hetero_m_p < 0.05:
            print(f"the data generated for {gene_symbol} for parameter {parameter_stable_id} is between centers for heterozygote males not very similar/comparable")
            
        if hetero_f_p < 0.05:
            print(f"the data generated for {gene_symbol} for parameter {parameter_stable_id} is between centers for heterozygote females not very similar/comparable")
        
        print("-----------------------------------------------------------------")
        
        
        #####make categorical plot
        
        plt.figure(figsize=(10, 12), dpi = 500)
    
        sns.set_theme()
    
        box = sns.catplot(data=df, x="phenotyping_center", y="knockout_wt_ratio", col="sex", row = "zygosity", kind = plot_type)
            #category plot
            
        #move the following title higher
        if len(set(df.zygosity)) == 1:
            box.fig.subplots_adjust(top=.85)
        else:
            box.fig.subplots_adjust(top=.92)
    
        box.fig.suptitle(f'Knockout-WT-Ratio for {parameter_name} ({parameter_stable_id}) of {gene_symbol}, week: min = {min_week}, max = {max_week}', size = 18, ha = 'center')
        #bug: ha/horizontalaligment= 'right' ----> 'left'
        
        
        for ax in box.axes.flat:
            ax.tick_params(axis='x', labelrotation = 90)
            ax.set_xlabel("Phenotyping Center", fontsize=10)
            ax.set_ylabel("Knockout-WT-Ratio", fontsize=10)
            
        if focus_zygosity == False:
            zygosity = "all_zygosities"
        
        if less_outliers == True:
            plt.savefig(f"{path}/{plot_type}plot_{gene_symbol}_{parameter_stable_id}_kwr_{zygosity}_week_min_{min_week}_max_{max_week}_less_outliers.png", dpi = 500, bbox_inches='tight')
            #whole path for plot saving important
            
            results.to_csv(f'{path}/{gene_symbol}_{parameter_stable_id}_ANOVA_{zygosity}_week_min_{min_week}_max_{max_week}_less_outliers.csv', index=False)
            
            if save_memory == False:
                globals()[f"{gene_symbol}_{parameter_stable_id}_ANOVA_{zygosity}_week_min_{min_week}_max_{max_week}_less_outliers"] = results
                # global dynamic variable gets results assigned
                    # no need to return a variable
                
            if len(all_post_hoc_p) != 0:
                all_post_hoc_p.to_csv(f'{path}/{gene_symbol}_{parameter_stable_id}_post_hoc_p_{post_hoc_method}_{zygosity}_week_min_{min_week}_max_{max_week}_less_outliers.csv', index=False)
                if save_memory == False:
                    globals()[f"{gene_symbol}_{parameter_stable_id}_post_hoc_p_{zygosity}_week_min_{min_week}_max_{max_week}_less_outliers"] = all_post_hoc_p
            
                
        if less_outliers == False:
            plt.savefig(f"{path}/{plot_type}plot_{gene_symbol}_{parameter_stable_id}_kwr_{zygosity}_week_min_{min_week}_max_{max_week}_outliers_incl.png", dpi = 500, bbox_inches='tight')
            #whole path for plot saving important
            
            results.to_csv(f'{path}/{gene_symbol}_{parameter_stable_id}_ANOVA_{zygosity}_week_min_{min_week}_max_{max_week}_outliers_incl.csv', index=False)
            
            if save_memory == False:
                globals()[f"{gene_symbol}_{parameter_stable_id}_ANOVA_{zygosity}_week_min_{min_week}_max_{max_week}_outliers_incl"] = results
                
            if len(all_post_hoc_p) != 0:
                all_post_hoc_p.to_csv(f'{path}/{gene_symbol}_{parameter_stable_id}_post_hoc_p_{post_hoc_method}_{zygosity}_week_min_{min_week}_max_{max_week}_outliers_incl.csv', index=False)
                if save_memory == False:
                    globals()[f"{gene_symbol}_{parameter_stable_id}_post_hoc_p_{zygosity}_week_min_{min_week}_max_{max_week}_outliers_incl"] = all_post_hoc_p
          
                
        
        if save_memory == True:
            plt.clf()
            #clear current figure
#%%
#call function
genes = ["Prkab1", "Dbn1", "Ap4e1", "Nxn"]
BW_ANOVA_center(genes, min_week = 14, max_week = 16, save_memory = True)
