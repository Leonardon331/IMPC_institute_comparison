# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 13:43:53 2023

@author: Leonard Borchardt
"""


import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
from PIL import Image

#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')
#%%

def overall_ANOVA_analysis(path: str = r'./parameters', focus_zygosity: bool =  True, less_outliers: bool = False):
    '''
    
    1. concatenates ANOVA data frames and post hoc test data frames
    
    2. creats a new data frame:
        counts relates how often in general and how often the phenotyping centers differ from the others
    
    plots:
     - pie chart with rel. parts of significant and non-significant ANOVA results
     - bar plot with rel. amount of significant ANOVAs and post hoc tests of each center

    import modules:
      - import os
      - import pandas as pd
      - import matplotlib.pyplot as plt
      - import seaborn as sns
      - import plotly.graph_objects as go
      - from PIL import Image
        
    Parameters
    ----------
     - path: str, cwd/parameters by default, optional
     
     these have to match the ANOVA_center.py/BW_ANOVA_center.py parameters:
     ----------------------------------------------------------------------
     - focus_zygosity: bool, default: True; if True: comparison focuses on zygosity with most data across centers, optional
     - less_outliers: bool, default: False, take only ANOVA results into count with or without strong (quartile+/-3xIQR) (default) outliers, optional


    Returns and Saves
    -------
     -None
     
    automatically assigns variables for:
     - data frame with absolute and relative values of significant ANOVAs and post hoc tests for each center
       
    '''
    
    
    save = r'./overall_analysis'
    # check whether directory already exists, otherwise create new folder
    if not os.path.exists(save):
        os.mkdir(save)
    
    save = r'./overall_analysis/ANOVAs'
    # check whether directory already exists, otherwise create new folder
    if not os.path.exists(save):
        os.mkdir(save)
    
    
    ### concatenate ANOVA and post hoc results ###
    print ("concatenate all ANOVA results")
    
    def add_ANOVA(less_outliers = less_outliers, focus_zygosity = focus_zygosity):
        #outer fct.
        
        contains = 'ANOVA' # substring that the filename should contain
        if less_outliers == False:
            outliers = "outliers_incl"
        else:
            outliers = "less_outliers"
        
        if focus_zygosity == False:
            zyg = "all_zygosities"
        else:
            zyg = "zygote"
        
        
        def add_ANOVA_list(path: str = r'./parameters', ANOVAs: list = []):
            #recursive fct.
            
            for item in os.listdir(path):
                item_path = os.path.join(path, item)
                if os.path.isfile(item_path) and contains in item and outliers in item and zyg in item:
                    df = pd.read_csv(item_path)
                    ANOVAs.append(df)
                    print(item, "appended")
                    
                if os.path.isdir(item_path):
                    add_ANOVA_list(item_path, ANOVAs)
                
            return ANOVAs
        
        ANOVAs = add_ANOVA_list()
        
        if focus_zygosity == True:
            zyg = "single_zygosities"
            
        save = r'./overall_analysis'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(save):
            os.mkdir(save)
        save = r'./overall_analysis/ANOVAs'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(save):
            os.mkdir(save)
        
        if len(ANOVAs) != 0:      
            ANOVAs = pd.concat(ANOVAs,ignore_index=True)
            
            ANOVAs = ANOVAs[~ANOVAs.p_Value_ANOVA.apply(lambda x: isinstance(x, str))]
                # Drop rows with string data in the 'p_Value_ANOVA' column
            
            ANOVAs = ANOVAs.sort_values(by=["p_Value_ANOVA"]).reset_index(drop=True)
            ANOVAs.to_csv(f'{save}/ANOVAs_{zyg}_{outliers}.csv', index=False)
            print("all ANOVA results concatenated")
            
        else:
            print('check if the right ANOVA data frames exist (e.g. for all zygosities with outliers:  "ANOVA_..._all_zygosities_outliers_incl")')
         
    add_ANOVA()
    
    
    print()
    print("concatenate all post hoc results")
    def add_post_hocs(less_outliers = less_outliers, focus_zygosity = focus_zygosity):
        #outer fct.
        
        contains = 'post_hoc_p'
        
        if less_outliers == False:
            outliers = "outliers_incl"
        else:
            outliers = "less_outliers"
        
        if focus_zygosity == False:
            zyg = "all_zygosities"
        else:
            zyg = "zygote"
        
        
        def add_post_hocs_list(path: str = r'./parameters', post_hocs: list = []):
            #recursive fct.   
            
            for item in os.listdir(path):
                item_path = os.path.join(path, item)
                if os.path.isfile(item_path) and contains in item and outliers in item and zyg in item:
                    df = pd.read_csv(item_path)
                    post_hocs.append(df)
                    print(item, "appended")
                    
                if os.path.isdir(item_path):
                    add_post_hocs_list(item_path, post_hocs)
        
            return post_hocs
        
        post_hocs = add_post_hocs_list()
        
        
        if focus_zygosity == True:
            zyg = "single_zygosities"
        
        save = r'./overall_analysis'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(save):
            os.mkdir(save)
        save = r'./overall_analysis/ANOVAs'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(save):
            os.mkdir(save)
        
        if len(post_hocs) != 0:      
            post_hocs = pd.concat(post_hocs,ignore_index=True)
            post_hocs.to_csv(f'{save}/post_hocs_{zyg}_{outliers}.csv', index=False)
            print("all post hoc results concatenated")
    
        else:
            print('check if the right ANOVA data frames exist (e.g. for all zygosities with outliers:  "post_hoc_p_..._all_zygosities_outliers_incl")')   
    
    add_post_hocs()
    
    
    if less_outliers == False:
        outliers = "outliers_incl"
    else:
        outliers = "less_outliers"
    
    if focus_zygosity == False:
        zyg = "all_zygosities"
    else:
        zyg = "single_zygosities"
    
    
    post_hocs = pd.read_csv(f'{save}/post_hocs_{zyg}_{outliers}.csv')
    
    
    print()
    print()
    print("make pie chart for proportion of overall ANOVA significance")
    ##### across all centers: relative ANOVA significance  #####
    
    ANOVAs = pd.read_csv(f'{save}/ANOVAs_{zyg}_{outliers}.csv')
    
    count_ANOVAs = len(ANOVAs.groupby(["gene_symbol","parameter_stable_id","zygosity","sex"])["p_Value_ANOVA"].mean())
    
    count_sig_ANOVAs = len(ANOVAs[ANOVAs.p_Value_ANOVA <= 0.05].groupby(["gene_symbol","parameter_stable_id","zygosity","sex"])["p_Value_ANOVA"].mean())
        
    count_non_sig_ANOVAs = count_ANOVAs - count_sig_ANOVAs
    
    #rel_sig_ANOVAs = count_sig_ANOVAs/count_ANOVAs
    #rel_non_sig_ANOVAs = count_non_sig_ANOVAs/count_ANOVAs
    
    
    #### make pie chart
    labels = ["Non-Significant ANOVAs", "Significant ANOVAs"]
    values = [count_non_sig_ANOVAs, count_sig_ANOVAs]
    
    fig = go.Figure(data=[go.Pie(labels=labels,
                                 values=values,
                                 texttemplate = "%{label}: <br> %{value} (%{percent})",
                                 textfont_size=20,
                                 marker=dict(colors=["mediumturquoise", "darkorange"]),
                                 pull=[0, 0.2])])
    
    fig.update_layout(showlegend = False,
                      title = dict(
                          text = "Proportion of significant ANOVAs",
                          font_size = 30,
                          automargin = True,
                          x = 0.5,
                          y = 0.9
                          ))
      
    fig.write_image(f"{save}/pie_chart_significant_ANOVAs_{zyg}_{outliers}.png", scale = 2, width=700, height=700)
    
    
    with Image.open(f"{save}/pie_chart_significant_ANOVAs_{zyg}_{outliers}.png") as pie:
        width, height = pie.size
        
        # Setting the points for cropped image
        left = width - width*0.95
        top = height - height*0.95
        right = width*0.95
        bottom = height*0.90
        
        pie = pie.crop((left, top, right, bottom))
        
        pie.save(f"{save}/pie_chart_significant_ANOVAs_{zyg}_{outliers}.png")
        pie.show()
        
        
        
    print()
    print()
    print("calculate relative ANOVA and post hoc significance of each center (center_statistics)")
    ######### center stats: relative significance of each center ########
    center_statistics = []
    
    def row_sum(row):
        count = 0
        for element in row:
            if isinstance(element, (int, float)) and element <= 0.05:
                count += 1
        return count


    for center in list(set(ANOVAs.phenotyping_center.sort_values(ignore_index=True))):
        
        ### ANOVA center stats ###
        
        df = ANOVAs[ANOVAs.phenotyping_center == center]
        num_ANOVAs_center = len(df)
        num_sig_ANOVAs_center = df.p_Value_ANOVA.to_frame().apply(lambda row: row_sum(row), axis=1).sum() ##########
        #print(num_sig_ANOVAs_center)
        rel_sig_ANOVAs_center = num_sig_ANOVAs_center/num_ANOVAs_center
    
        

        ### post hoc center stats ###
    
        center_ph = post_hocs[post_hocs.phenotyping_center==center]

        # Count values smaller/equal 0.05 for the center
        absl_significant = center_ph.apply(lambda row: row_sum(row), axis=1).sum()
    
        ###all possible post hoc tests with the center:
        all_tests = 0
        for parameter_stable_id in set(ANOVAs.parameter_stable_id):
            for zygosity in set(ANOVAs[ANOVAs.parameter_stable_id == parameter_stable_id].zygosity):
                for sex in set(ANOVAs[(ANOVAs.parameter_stable_id == parameter_stable_id) & (ANOVAs.zygosity == zygosity)].sex):
                    if center in set(ANOVAs[(ANOVAs.parameter_stable_id == parameter_stable_id) & (ANOVAs.zygosity == zygosity) & (ANOVAs.sex == sex)].phenotyping_center):
                        all_tests += len(ANOVAs[(ANOVAs.parameter_stable_id == parameter_stable_id) & (ANOVAs.zygosity == zygosity) & (ANOVAs.sex == sex)])-1
                                #-1 bcs. no self test possible
         
        
        rel_significant = absl_significant/all_tests
        
        center_stats = pd.DataFrame({
            "phenotyping_center": center,
            "ANOVAs": [num_ANOVAs_center],
            "significant_ANOVAs": [num_sig_ANOVAs_center],
            "rel_significant ANOVAs": [rel_sig_ANOVAs_center],
            "all_possible_post_hocs": [all_tests],
            "significant_post_hocs": [absl_significant],
            "rel_significant_post_hocs": [rel_significant]
            })
        
        center_statistics.append(center_stats)
        
    center_statistics = pd.concat(center_statistics, ignore_index=True)
    
    center_statistics = center_statistics.sort_values(by=['phenotyping_center']).reset_index(drop=True)

    center_statistics.to_csv(f'{save}/center_statistics_{zyg}_{outliers}.csv', index=False)
    
    globals()["center_statistics_{zyg}_{outliers}"] = center_statistics
    print("center statistics successfully calculated")
    
    
    print()
    print()
    print("make bar plot for comparing the significance of each center")
    ### barplot ###    
    plt.figure(figsize=(8, 8), dpi = 350)

    sns.set_theme()
    
    ########## double boxplot, with significant ANOVAS
    
    ##bar = center_statistics.plot(x="phenotyping_center", y=["rel_significant ANOVAs", "rel_significant_post_hocs"], kind="bar")
    #bar = sns.barplot(data=center_statistics, x="phenotyping_center", y=["rel_significant ANOVAs","rel_significant_post_hocs"])
                                                                        #nur 1 y bei sns möglich

    ##bar.figure.subplots_adjust(top=.87)


    ##bar.figure.suptitle('Relative number of significant ANOVAs and post hoc tests per center', size = 18, ha = 'center')
    #bug: ha/horizontalaligment= 'right' ----> 'left'
    
    
    ##bar.tick_params(axis='x', labelrotation = 45)
    ##bar.set_xlabel("Phenotyping Center", fontsize=14)
    ##bar.set_ylabel("Relative Significance", fontsize=14)
    ##plt.legend(["Significant ANOVAs", "Significant Post Hocs"],bbox_to_anchor = (1.5, 0.6), loc='center right')
    
    #plt.savefig(f"{save}/bar_plot_rel_ANOVA_and_post_hoc_center_significance_{zyg}_{outliers}.png", dpi = 500, bbox_inches='tight')
    
    ####################
    
   
    bar = sns.barplot(data=center_statistics, x="phenotyping_center", y="rel_significant_post_hocs")
                                                                        
    bar.figure.subplots_adjust(top=.92)

    bar.figure.suptitle('Relative number significant post hoc tests per center', size = 18, ha = 'center')
    #bug: ha/horizontalaligment= 'right' ----> 'left'
    
    bar.tick_params(axis='x', labelrotation = 90)
    bar.set_xlabel("Phenotyping Center", fontsize=14)
    bar.set_ylabel("Relative Significance", fontsize=14)
   
   
    #for i in bar.containers:
        #bar.bar_label(i,)
        # rel. value for center above bar
    
    plt.savefig(f"{save}/bar_plot_rel_post_hoc_center_significance_{zyg}_{outliers}.png", dpi = 500, bbox_inches='tight')
    #whole path for plot saving important
    
#%%
#call function
overall_ANOVA_analysis()