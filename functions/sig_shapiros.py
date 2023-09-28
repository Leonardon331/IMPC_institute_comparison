# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 17:41:10 2023

@author: Leonard Borchardt
"""

import os
import warnings
import pandas as pd
import scipy.stats as stats
import plotly.graph_objects as go
from PIL import Image

#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')
#%%

def sig_shapiros():
    '''
    creates data frame with knockout-wt-ratio group (gene/sex/zygosity/parameter/institute) shapiro test
    significance proportion (of all tests)
    creates pie chart
    
    import modules:
        import os
        import warnings
        import pandas as pd
        import scipy.stats as stats
        import plotly.graph_objects as go
        from PIL import Image
    
    run: concat_knockout_wt_ratios.py first
    
    Parameters
    ----------
    None

    Returns
    -------
    None
    
    creates variable with data frame results
    saves the dataframe and the pie chart under cwd/overall_analysis/kwr_normality

    '''
    
    # ignore warnings
    warnings.filterwarnings("ignore")
    
    save = r'./overall_analysis'
    # check whether directory already exists, otherwise create new folder
    if not os.path.exists(save):
        os.mkdir(save)
    
    save = r'./overall_analysis/kwr_normality'
    # check whether directory already exists, otherwise create new folder
    if not os.path.exists(save):
        os.mkdir(save)
    
    read = 'overall_analysis/knockout_wt_ratio/cleaned_knockout_wt_ratios.csv'
    df = pd.read_csv(read)
    
    groups = df.groupby(["gene_symbol","parameter_stable_id","phenotyping_center","sex","zygosity"])
    
    print("perform group Shapiro Tests..")
    
    shapiros = []
    for name, group in groups:
        if len(group) < 3:
            continue
        W,p = stats.shapiro(group.knockout_wt_ratio)

        shapiro = pd.DataFrame(dict(
            gene_symbol = name[0],
            procedure_stable_id = group.procedure_stable_id.unique(),
            procedure_name = group.procedure_name.unique()[0],  # 0 bcs. of occuring multiple procedure names
            parameter_stable_id = name[1],
            parameter_name = group.parameter_name.unique()[0],
            phenotyping_center = name[2],
            sex = name[3],
            zygosity = name[4],
            shapiro_W = W,
            shapiro_p = p,
            sample_size = len(group)
            ))
        
        
        shapiros.append(shapiro)
        
        
    shapiros = pd.concat(shapiros,ignore_index=True)
    print("all Shapiro Tests performed")
    
    print("calculate sginificance proportion and make pie chart..")
    shapiro_counts = pd.DataFrame([len(shapiros)],columns = ["absl_all"])
    shapiro_counts["absl_sig"] = [len(shapiros[shapiros.shapiro_p<=0.05])]
    shapiro_counts["rel_sig"] = shapiro_counts.absl_sig.to_numpy()/shapiro_counts.absl_all.to_numpy()
    
    
    # make pie chart
    labels = ["Non-Significant Shapiros", "Significant Shapiros"]
    values = [shapiro_counts.absl_all.to_numpy()[0]-shapiro_counts.absl_sig.to_numpy()[0], shapiro_counts.absl_sig.to_numpy()[0]]
    
    fig = go.Figure(data=[go.Pie(labels=labels,
                                 values=values,
                                 texttemplate = "%{label}: <br> %{value} (%{percent})",
                                 textfont_size=20,
                                 marker=dict(colors=["mediumturquoise", "darkorange"]),
                                 pull=[0, 0.2])])
    
    fig.update_layout(showlegend = False,
                      title = dict(
                          text = "Proportion of Significant Shapiro Tests",
                          font_size = 30,
                          automargin = True,
                          x = 0.5,
                          y = 0.9
                          ))
      
    fig.write_image(f"{save}/pie_chart_significant_Shapiros.png", scale = 2, width=700, height=700)
    
    
    with Image.open(f"{save}/pie_chart_significant_Shapiros.png") as pie:
        width, height = pie.size
        
        # Setting the points for cropped image
        left = width - width*0.95
        top = height - height*0.95
        right = width*0.95
        bottom = height*0.90
        
        pie = pie.crop((left, top, right, bottom))
        
        pie.save(f"{save}/pie_chart_significant_Shapiros.png")
        pie.show()
    
    
    
    shapiros = shapiros.sort_values(by=["shapiro_p"]).reset_index(drop=True)
    shapiros.to_csv(f"{save}/shapiros.csv", index = False)
    shapiro_counts.to_csv(f"{save}/shapiro_counts.csv", index = False)
    globals()["shapiro_counts"] = shapiro_counts
    print("All done.")
#%%
sig_shapiros()