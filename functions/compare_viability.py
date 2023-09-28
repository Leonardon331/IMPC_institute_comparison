# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 15:21:16 2023

@author: Leonard Borchardt
"""

import os
import warnings
import pandas as pd
import numpy as np
import altair as alt


#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')
#%%

def compare_viability(gene_symbols: list):
    '''
    creates a concatenated data frame with homozygote viability for every gene, center and sex
    creates a bar with the viability outcome proportions for all genes
    
    import modules:
        import os
        import warnings
        import pandas as pd
        import numpy as np
        import altair as alt
    
    dowload genes first with download_genes.py
    
    Parameters
    ----------
    gene_symbols: list of strings, gene symbols for genes downloaded

    Returns
    -------
    None
    
    automatically assigns variable for gene data frame with viability-center results
    saves the data frame , too
    '''
    
    #ignore warnings
    warnings.filterwarnings("ignore")
    
    print()
    print()
    
    viabilities = []
    for gene_symbol in gene_symbols:
        print(gene_symbol)
        read = f'./genes/{gene_symbol}/{gene_symbol}.csv'   ##read in the gene data
        df = pd.read_csv(read)
    
        ########
        #create folder if it doesn't already exist
        save = r'./overall_analysis'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(save):
            os.mkdir(save)
        
        #create subfolder if it doesn't already exist
        save = f'{save}/viability'
        # check whether directory already exists, otherwise create new folder
        if not os.path.exists(save):
            os.mkdir(save)
        
        #######
        
        # only viability data
        df = df[df.procedure_name == "Viability Primary Screen"]
        
        #drop nan values
        df = df.dropna(subset=['data_point'])
        
        centers = list(set(df.phenotyping_center))
        centers.sort()
        
        #choosing the ps_ids according to center measurement pattern (observed with gene df filtering in excel)
        for center in centers:
            df_center = df[df.phenotyping_center == center]
            colonies = df_center.groupby("colony_id") # group after colony_id to remove double measurements for one colony
            sexes_df = []
            for sex in ["male", "female"]:
                if sex == "male":
                    #"IMPC_VIA_007_001": num. male wts
                    #"IMPC_VIA_008_001": num. male het
                    #"IMPC_VIA_009_001": num. male homo
                    #"IMPC_VIA_049_001": num. male wts
                    #"IMPC_VIA_051_001": num. male het
                    #"IMPC_VIA_053_001": num. male homo
                    
                    if df_center.parameter_stable_id.isin(["IMPC_VIA_007_001","IMPC_VIA_008_001","IMPC_VIA_009_001","IMPC_VIA_049_001","IMPC_VIA_051_001","IMPC_VIA_053_001"]).any():
                        wt = []
                        het = []
                        hom = []
                        hemi = []
                        for name, group in colonies:
                            #center colony male means
                            col_wt_m = group[group.parameter_stable_id.isin(["IMPC_VIA_007_001", "IMPC_VIA_049_001"])].data_point.mean()
                            wt.append(col_wt_m)
                            col_het_m = group[group.parameter_stable_id.isin(["IMPC_VIA_008_001", "IMPC_VIA_051_001"])].data_point.mean()
                            het.append(col_het_m)
                            col_hom_m = group[group.parameter_stable_id.isin(["IMPC_VIA_009_001", "IMPC_VIA_053_001"])].data_point.mean()
                            hom.append(col_hom_m)

                            if group.parameter_name.isin(["Total of hemizygous males" ,"Total male hemizygous"]).any():      #hemizygotes are always males
                                col_hemi_m = group[group.parameter_name.isin(["Total of hemizygous males", "Total male hemizygous"])].data_point.mean()
                                hemi.append(col_hemi_m)
                            else:
                                hemi.append(0)
                        
                    #sum of all colony means
                    wt = np.array(wt).sum()
                    het = np.array(het).sum()
                    hom = np.array(hom).sum()
                    hemi = np.array(hemi).sum()
                    wt = np.array(wt).sum()
                    rel_homo = hom/(hom+het+hemi+wt)
                    absl_pups = het + hom + wt + hemi
                    absl_homo = rel_homo*absl_pups
                   
                    if absl_pups == 0:
                        continue  #next sex
                    
                if sex == "female":
                    #"IMPC_VIA_011_001": num. female wts
                    #"IMPC_VIA_012_001": num. female het
                    #"IMPC_VIA_013_001": num. female homo
                    #"IMPC_VIA_050_001": num. female wts
                    #"IMPC_VIA_052_001": num. female het
                    #"IMPC_VIA_054_001": num. female homo
                    
                    if df_center.parameter_stable_id.isin(["IMPC_VIA_011_001","IMPC_VIA_012_001","IMPC_VIA_013_001", "IMPC_VIA_050_001","IMPC_VIA_052_001","IMPC_VIA_054_001"]).any():
                        wt = []
                        het = []
                        hom = []
                        for name, group in colonies:
                            #center colony females het/hom
                            col_wt_f = group[group.parameter_stable_id.isin(["IMPC_VIA_011_001", "IMPC_VIA_050_001"])].data_point.mean()
                            wt.append(col_wt_f)
                            col_het_m = group[group.parameter_stable_id.isin(["IMPC_VIA_012_001", "IMPC_VIA_052_001"])].data_point.mean()
                            het.append(col_het_m)
                            col_hom_m = group[group.parameter_stable_id.isin(["IMPC_VIA_013_001", "IMPC_VIA_054_001"])].data_point.mean()
                            hom.append(col_hom_m)
                    
                    het = np.array(het).sum()
                    hom = np.array(hom).sum()
                    wt = np.array(wt).sum()
                    rel_homo = hom/(hom+het+wt)
                    absl_pups = het + hom + wt
                    absl_homo = rel_homo*absl_pups
                    
                    if absl_pups == 0:
                        continue  #next sex
            
                    
                if rel_homo >= 0.12:
                    category = "Homozygous - Viable"
                elif rel_homo < 0.12 and rel_homo > 0:
                    category = "Homozygous - Subviable"
                elif rel_homo == 0:
                    category = "Homozygous - Lethal"
                else:
                    category = float("nan")
                
                viability = pd.DataFrame(dict(
                    gene_symbol = [gene_symbol],
                    sex = [sex],
                    phenotyping_center = [center],
                    absl_pups = [absl_pups],
                    absl_homo = [absl_homo],
                    rel_homozygous = [rel_homo],
                    category = [category]
                    ))
                sexes_df.append(viability)
                viabilities.append(viability)
            
            sexes_df = pd.concat(sexes_df, ignore_index=True) #temp data frame 
            
            
            sex = "both"
            #append "both" sex values
            #absl homozygotes for all sexes
            hom = [sexes_df.at[0,"rel_homozygous"]*sexes_df.at[0,"absl_pups"], #male
                   sexes_df.at[1,"rel_homozygous"]*sexes_df.at[1,"absl_pups"]] #female
            #all offspring
            absl_pups = [sexes_df.at[0,"absl_pups"], sexes_df.at[1,"absl_pups"]] #row,column ---> access one value in df
            
            hom = np.array(hom).sum()   #sum of colonies ---> all hom animals at center
            absl_pups = np.array(absl_pups).sum()
            rel_homo = hom/absl_pups
            absl_homo = rel_homo*absl_pups
            
            if absl_pups == 0:
                continue  #next center
                                                 
            if rel_homo >= 0.12:
                category = "Homozygous - Viable"
            elif rel_homo < 0.12 and rel_homo > 0:
                category = "Homozygous - Subviable"
            elif rel_homo == 0:
                category = "Homozygous - Lethal"
            else:
                category = float("nan")
            
            viability = pd.DataFrame(dict(
                gene_symbol = [gene_symbol],
                sex = [sex],
                phenotyping_center = [center],
                absl_pups = [absl_pups],
                absl_homo = [absl_homo],
                rel_homozygous = [rel_homo],
                category = [category]
                ))
            viabilities.append(viability)
            
        print(f"center viabilities for {gene_symbol} calculated")
        
    print()
    
    viabilities = pd.concat(viabilities, ignore_index=True)
    viabilities = viabilities.sort_values(by=["gene_symbol","sex","phenotyping_center"], ascending=[True,False,True]).reset_index(drop=True)
    viabilities.to_csv(f"{save}/viability.csv", index=False) 
    globals()["viabilities"] = viabilities
    print("viabilities combined")
    
    
    ## count viability categories occurence for eache gene and sex
    viabilities_count = viabilities.groupby(["gene_symbol","sex","category"]).size().reset_index(name='count')

    # add relative values, category occurance regarding all gene/sex measurements/classifications
    genes_count = viabilities_count.groupby(["gene_symbol","sex"])
    rel_cat_counts = []
    for name,group in genes_count:
        all_cat_count = group["count"].sum()    #gene/sex classification count
                                        #np.array with single element
        rel_cat_count = list(np.array(viabilities_count[(viabilities_count.gene_symbol == name[0])&(viabilities_count.sex == name[1])]["count"])/all_cat_count)
                                                                                                        #access category count series for gene and sex of all_cat_count
        rel_cat_counts += rel_cat_count # fill list/column
    viabilities_count.insert(4, "rel_count", rel_cat_counts)
    
    # implement list of centers for categories
    centers = viabilities.groupby(["gene_symbol","sex","category"])["phenotyping_center"].apply(list).reset_index(name='phenotyping_centers')
    viabilities_count.insert(5,"phenotyping_centers",centers["phenotyping_centers"])
    
    viabilities_count.category=pd.Categorical(viabilities_count.category,categories=["Homozygous - Viable", "Homozygous - Subviable", "Homozygous - Lethal"])  #force this specific order
    viabilities_count = viabilities_count.sort_values(by=["category"]).reset_index(drop=True) # sort in the specific order
    viabilities_count = viabilities_count.sort_values(by=["gene_symbol","sex"], ascending=[True,False]).reset_index(drop=True)
    
    viabilities_count.to_csv(f"{save}/viability_count.csv", index=False)
    globals()["viability_count"] = viabilities_count
    print("viabilities counted")



    ### stacked viability bar chart with bars for each sex and gene ###
    print()
    print("make bar chart with proportions of viability categories for all genes")
    
    #following this approach:
     #https://stackoverflow.com/questions/22787209/how-to-have-clusters-of-stacked-bars

    
    viable = viabilities_count[viabilities_count.category == "Homozygous - Viable"][["gene_symbol", "sex","rel_count"]].reset_index(drop=True)
    subviable = viabilities_count[viabilities_count.category == "Homozygous - Subviable"][["gene_symbol", "sex","rel_count"]].reset_index(drop=True)
    lethal = viabilities_count[viabilities_count.category == "Homozygous - Lethal"][["gene_symbol", "sex","rel_count"]].reset_index(drop=True)
    
    # append zero rel_counts
    def prep_df(df, category, gene_symbols = gene_symbols):
        #append missing gene/sex rows to the end of df
        gene_symbols.sort()
        sexes = ["male","female","both"]
        rows = []
        for gene_symbol in gene_symbols:
            for sex in sexes:
                if df[(df.gene_symbol == gene_symbol)&(df.sex == sex)].empty:
                    row = [gene_symbol,sex,0] #zero rel_count row to append at the end of df
                    rows.append(row)
        rows = pd.DataFrame(rows, columns = ["gene_symbol","sex","rel_count"]) # make list of lists to df           
        #append missing rows to df   
        df = pd.concat([df,rows],ignore_index=True) # append missing rows to df
        df = df.sort_values(by=["gene_symbol","sex"], ascending=[True,False]).reset_index(drop=True)
        
        df["category"] = category
        return df
    
    # utilize function to append rows
    viable = prep_df(viable, 'Viable')
    subviable = prep_df(subviable, 'Subviable')
    lethal = prep_df(lethal, 'Lethal')
    df = pd.concat([viable, subviable, lethal])  #concat data frames
    
    
    #Plot data with Altair
    chart = alt.Chart(df).mark_bar().encode(

        # tell Altair which field to group columns on
        x=alt.X('sex:N', title=None,
                sort=alt.EncodingSortField(field="sex",order='descending')),
    
        # tell Altair which field to use as Y values and how to calculate
        y=alt.Y('sum(rel_count):Q',
            axis=alt.Axis(
                grid=False,
                title="Proportion")),
    
        # tell Altair which field to use to use as the set of columns to be  represented in each group
        column=alt.Column('gene_symbol:N', title=None),
    
        # tell Altair which field to use for color segmentation 
        color=alt.Color('category:N', title = None,
                scale=alt.Scale(
                    # make it look pretty with an enjoyable color pallet
                    range=['#ff6f69','#ffcc5c','#96ceb4'],
                ),
            ))\
        .configure_view(
            # remove grid lines around column clusters
            strokeOpacity=0    
        ).properties(
    title="Homozygote Viability Proportion"
    )
     
    chart.save(f'{save}/bar_chart_viability.png', ppi=500, scale_factor=1.5)
#%%
#run function
compare_viability(gene_symbols = ["Prkab1","Dbn1","Ap4e1","Nxn"])