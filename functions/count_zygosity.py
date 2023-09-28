# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 09:38:04 2023

@author: Leonard Borchardt
"""

import os
import warnings
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')
#%%
def count_zygosity():
    '''
    counts for each gene the zygosities of the mice tested and plots the proportions
    
    import modules:
        import os
        import warnings
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns
        
    needs concatenated gene data frame with all columns,
    run:
        "download_genes.py"
        "concat_genes.py"
    
    Parameters
    ----------
    None

    Returns
    -------
    None

    automatically assigns variable for data frame with gene zygosity counts
    saves the data frame under cwd/overall_analysis_zygosity/zygosity_counts.csv
    saves also the proportions and a stacked bar plot with the proportions
    '''
    
    # ignore warnings
    warnings.filterwarnings("ignore")
    
    ########
    #create folder if it doesn't already exist
    save = r'./overall_analysis'
    # check whether directory already exists, otherwise create new folder
    if not os.path.exists(save):
        os.mkdir(save)
    
    #create subfolder for saving if it doesn't already exist
    save = f'{save}/zygosity'
    # check whether directory already exists, otherwise create new folder
    if not os.path.exists(save):
        os.mkdir(save)
    #########
    
    print("count zygosities..")
    
    df = pd.read_csv(r'./genes_all_columns/all_genes_all_columns.csv')
    
    #exclude viability offspring
    df = df[~df.procedure_stable_id.str.contains("VIA")]
    
    gr = df.groupby(["gene_symbol","zygosity"])["specimen_id"].nunique().reset_index(name="count")
    #.nunique(): counts unique values, leaves nan values out per default
    
    gene_symbols = list(set(gr.gene_symbol))
    gene_symbols.sort()

    rows = []
    for gene_symbol in gene_symbols:
        row = []
        row.append(gene_symbol)
        for zygosity in ["homozygote","heterozygote","hemizygote"]:
            z_count = gr[(gr.gene_symbol == gene_symbol)&(gr.zygosity == zygosity)]["count"]
            if z_count.empty:
                row.append(0)
            else:
              row.append(z_count.to_numpy()[0])  
        rows.append(row)
        
    zygosity_count = pd.DataFrame(rows, columns = ["gene_symbol","homozygote","heterozygote","hemizygote"])
    zygosity_count.set_index('gene_symbol', inplace = True)
    print("zygosities counted")
    zygosity_proportions = zygosity_count.apply(lambda x: x/sum(x), axis=1) # proportions in each row
    
    
    #stacked zygosity bar chart for each gene
    print("make stacked bar chart..")
    #adjust the subplot layout
    subplot_height = 8  # adjust this value to control the height
    subplot_width = 8   # adjust this value to control the width
    
    subplot_rows = 1
    subplot_columns = 1
    
    total_fig_height = subplot_height * subplot_rows
    total_figure_width = subplot_width * subplot_columns
    
    # Create subplots with adjusted figsize
    fig, axs = plt.subplots(subplot_rows, subplot_columns, figsize=(total_figure_width, total_fig_height), dpi=500)
   
    sns.set_theme()              
        
    zygosity_proportions.plot(kind='bar', stacked=True, ax = axs, legend=False)
    axs.set_xticklabels(gene_symbols)
    axs.tick_params(axis='x', labelrotation = 0)
    axs.set_xlabel("Silenced Gene", fontsize=14)
    axs.set_ylabel("Proportion", fontsize=14)

    #legend
    plt.legend(bbox_to_anchor=(1.02,1), loc='upper left', title = "Zygosity")
    
    
    fig.subplots_adjust(top=.85)
    fig.suptitle("Proportions of Mouse Zygosities for each Knockout", size = 18, ha = 'center')
    plt.tight_layout()
    
    plt.savefig(f'{save}/barplot_zygosity_proportions.png', dpi = 500, bbox_inches='tight')
    print("stacked bar chart created")
    
    zygosity_count.to_csv(f"{save}/zygosity_count.csv",index = True)
    zygosity_proportions.to_csv(f"{save}/zygosity_proportions.csv",index = True)
    globals()["zygosity_count"] = zygosity_count
    globals()["zygosity_proportions"] = zygosity_proportions
    print("all completed")

#%%
count_zygosity()