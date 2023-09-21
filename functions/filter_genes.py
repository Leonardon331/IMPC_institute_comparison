# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 17:52:08 2023

@author: Leonard Borchardt
"""
import pandas as pd
import os

#set wd
os.chdir(r"C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis")
#%%
def filter_genes(gene_symbols: list, parameter_stable_ids: list):
    
    '''
    filters gene data frames for parameter measurements
    
    Parameters:
        gene_symbols: list of strings, list of gene_symbols that were downloaded and saved as csv files via function download_genes
        parameter_stable_ids: list of strings, list of parameter measurement ids for filtering the gene data frames
        
    Returns:
        None
        automatically assigns variables to filtered gene data frames
    '''
    
    for gene_symbol in gene_symbols:
        path = f'./genes/{gene_symbol}/{gene_symbol}.csv'
        if os.path.exists(path):
            df = pd.read_csv(path) 
        else:
            print("download genes first")
            continue
        
        for parameter_stable_id in parameter_stable_ids:
            
            ######## write to:
            #create folder if it doesn't already exist
            path = f'./parameters/{parameter_stable_id}'
            # check whether directory already exists, otherwise create new folder
            if not os.path.exists(path):
                os.mkdir(path)
            
            #create subfolder for saving if it doesn't already exist
            path = f'.{path}/{gene_symbol}'
            # check whether directory already exists, otherwise create new folder
            if not os.path.exists(path):
                os.mkdir(path)
            #########
            
            df = df[df.parameter_stable_id == parameter_stable_id]
    
            df.to_csv(f'{path}/{gene_symbol}_{parameter_stable_id}.csv', index=False)
            globals()[f"{gene_symbol}_{parameter_stable_id}"] = df
                        #dynamic variable name = filtered gene data frame, to have the data frames in the variable explorer
#%%
filter_genes(["Prkab1", "Dbn1", "Ap4e1", "Nxn"], ["IMPC_ECG_028_001", "IMPC_XRY_001_001", "IMPC_EYE_007_001", "IMPC_CBC_025_001", "IMPC_CBC_011_001"])


















#%%
#Prkab1

Prkab1 = pd.read_csv("Prkab1.csv")

#body weight (curve)
Prkab1_IMPC_BWT_008_001 = Prkab1[Prkab1.parameter_stable_id == 'IMPC_BWT_008_001']
    #Plot weight curve for gene/sex and each center and compare data

#fertility
Prkab1_IMPC_FER_019_001 = Prkab1[Prkab1.parameter_stable_id == 'IMPC_FER_019_001']
Prkab1_IMPC_FER_001_001 = Prkab1[Prkab1.parameter_stable_id == 'IMPC_FER_001_001']

#viability
Prkab1_IMPC_VIA_002  = Prkab1[Prkab1.parameter_stable_id == 'IMPC_VIA_002']
    #needs to be calculated with formula (excel?)
Prkab1_IMPC_VIA_001  = Prkab1[Prkab1.parameter_stable_id == 'IMPC_VIA_001']

#ECG
Prkab1_IMPC_ECG_028_001 = Prkab1[Prkab1.parameter_stable_id == 'IMPC_ECG_028_001']
    #QT
Prkab1_IMPC_ECG_007_001 = Prkab1[Prkab1.parameter_stable_id == 'IMPC_ECG_007_001']
    #QRS

#IPGTT (Intraperitoneal glucose tolerance test)
Prkab1_IMPC_IPG_012_001 = Prkab1[Prkab1.parameter_stable_id == 'IMPC_IPG_012_001']
    #Area under glucose response curve (AUC)

#X-ray
Prkab1_IMPC_XRY_001_001 = Prkab1[Prkab1.parameter_stable_id == 'IMPC_XRY_001_001']
    #skull shape

#Eye
Prkab1_IMPC_EYE_007_001 = Prkab1[Prkab1.parameter_stable_id == 'IMPC_EYE_007_001']
    #cornea
Prkab1_IMPC_EYE_016_001 = Prkab1[Prkab1.parameter_stable_id == 'IMPC_EYE_0016_001']
    #lens
Prkab1_IMPC_EYE_017_001 = Prkab1[Prkab1.parameter_stable_id == 'IMPC_EYE_017_001']
    #lens opacity

#heart weight
Prkab1_IMPC_HWT_008_001 = Prkab1[Prkab1.parameter_stable_id == 'IMPC_HWT_008_001']
    #heart weight
Prkab1_IMPC_HWT_002_001 = Prkab1[Prkab1.parameter_stable_id == 'IMPC_HWT_002_001']
    #tibia length




#clinical chemistry
Prkab1_IMPC_CBC_025_001 = Prkab1[Prkab1.parameter_stable_id == 'IMPC_CBC_025_001']
    #LDL cholesterol

##iron
Prkab1_Fe = Prkab1[Prkab1.parameter_stable_id == 'IMPC_CBC_011_001']
##glucose
Prkab1_Glc = Prkab1[Prkab1.parameter_stable_id == 'IMPC_CBC_018_001']
##triglycerides
Prkab1_Tgc = Prkab1[Prkab1.parameter_stable_id == 'IMPC_CBC_017_001']
##potassium
Prkab1_K = Prkab1[Prkab1.parameter_stable_id == 'IMPC_CBC_002_001']
##sodium
Prkab1_Na = Prkab1[Prkab1.parameter_stable_id == 'IMPC_CBC_001_001']
##calcium
Prkab1_Ca = Prkab1[Prkab1.parameter_stable_id == 'IMPC_CBC_009_001']
##urea
Prkab1_Urea = Prkab1[Prkab1.parameter_stable_id == 'IMPC_CBC_011_001']
##total protein
Prkab1_Prot = Prkab1[Prkab1.parameter_stable_id == 'IMPC_CBC_006_001']
#%%
#Dbn1
Dbn1 = pd.read_csv("Dbn1.csv")

#body weight (curve)
Prkab1_IMPC_BWT_008_001 = Dbn1[Dbn1.parameter_stable_id == 'IMPC_BWT_008_001']
    #Plot weight curve for gene/sex and each center and compare data

#fertility
Dbn1_IMPC_FER_019_001 = Dbn1[Dbn1.parameter_stable_id == 'IMPC_FER_019_001']
Dbn1_IMPC_FER_001_001 = Dbn1[Dbn1.parameter_stable_id == 'IMPC_FER_001_001']

#viability
Dbn1_IMPC_VIA_002  = Dbn1[Dbn1.parameter_stable_id == 'IMPC_VIA_002']
    #needs to be calculated with formula (excel?)
Dbn1_IMPC_VIA_001  = Dbn1[Dbn1.parameter_stable_id == 'IMPC_VIA_001']

#ECG
Dbn1_IMPC_ECG_028_001 = Dbn1[Dbn1.parameter_stable_id == 'IMPC_ECG_028_001']
    #QT
Dbn1_IMPC_ECG_007_001 = Dbn1[Dbn1.parameter_stable_id == 'IMPC_ECG_007_001']
    #QRS

#IPGTT (Intraperitoneal glucose tolerance test)
Dbn1_IMPC_IPG_012_001 = Dbn1[Dbn1.parameter_stable_id == 'IMPC_IPG_012_001']
    #Area under glucose response curve (AUC)

#X-ray
Dbn1_IMPC_XRY_001_001 = Dbn1[Dbn1.parameter_stable_id == 'IMPC_XRY_001_001']
    #skull shape

#Eye
Dbn1_IMPC_EYE_007_001 = Dbn1[Dbn1.parameter_stable_id == 'IMPC_EYE_007_001']
    #cornea
Dbn1_IMPC_EYE_016_001 = Dbn1[Dbn1.parameter_stable_id == 'IMPC_EYE_0016_001']
    #lens
Dbn1_IMPC_EYE_017_001 = Dbn1[Dbn1.parameter_stable_id == 'IMPC_EYE_017_001']
    #lens opacity

#heart weight
Dbn1_IMPC_HWT_008_001 = Dbn1[Dbn1.parameter_stable_id == 'IMPC_HWT_008_001']
    #heart weight
Dbn1_IMPC_HWT_002_001 = Dbn1[Dbn1.parameter_stable_id == 'IMPC_HWT_002_001']
    #tibia length




#clinical chemistry
Dbn1_IMPC_CBC_025_001 = Dbn1[Dbn1.parameter_stable_id == 'IMPC_CBC_025_001']
    #LDL cholesterol

##iron
Dbn1_Fe = Dbn1[Dbn1.parameter_stable_id == 'IMPC_CBC_011_001']
##glucose
Dbn1_Glc = Dbn1[Dbn1.parameter_stable_id == 'IMPC_CBC_018_001']
##triglycerides
Dbn1_Tgc = Dbn1[Dbn1.parameter_stable_id == 'IMPC_CBC_017_001']
##potassium
Dbn1_K = Dbn1[Dbn1.parameter_stable_id == 'IMPC_CBC_002_001']
##sodium
Dbn1_Na = Dbn1[Dbn1.parameter_stable_id == 'IMPC_CBC_001_001']
##calcium
Dbn1_Ca = Dbn1[Dbn1.parameter_stable_id == 'IMPC_CBC_009_001']
##urea
Dbn1_Urea = Dbn1[Dbn1.parameter_stable_id == 'IMPC_CBC_011_001']
##total protein
Dbn1_Prot = Dbn1[Dbn1.parameter_stable_id == 'IMPC_CBC_006_001']

#%%
#Ap4e1

Ap4e1 = pd.read_csv("Ap4e1.csv")

#body weight (curve)
Ap4e1_IMPC_BWT_008_001 = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_BWT_008_001']
    #Plot weight curve for gene/sex and each center and compare data

#fertility
Ap4e1_IMPC_FER_019_001 = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_FER_019_001']
Ap4e1_IMPC_FER_001_001 = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_FER_001_001']

#viability
Ap4e1_IMPC_VIA_002  = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_VIA_002']
    #needs to be calculated with formula (excel?)
Ap4e1_IMPC_VIA_001  = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_VIA_001']

#ECG
Ap4e1_IMPC_ECG_028_001 = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_ECG_028_001']
    #QT
Ap4e1_IMPC_ECG_007_001 = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_ECG_007_001']
    #QRS

#IPGTT (Intraperitoneal glucose tolerance test)
Ap4e1_IMPC_IPG_012_001 = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_IPG_012_001']
    #Area under glucose response curve (AUC)

#X-ray
Ap4e1_IMPC_XRY_001_001 = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_XRY_001_001']
    #skull shape

#Eye
Ap4e1_IMPC_EYE_007_001 = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_EYE_007_001']
    #cornea
Ap4e1_IMPC_EYE_016_001 = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_EYE_0016_001']
    #lens
Ap4e1_IMPC_EYE_017_001 = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_EYE_017_001']
    #lens opacity

#heart weight
Ap4e1_IMPC_HWT_008_001 = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_HWT_008_001']
    #heart weight
Ap4e1_IMPC_HWT_002_001 = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_HWT_002_001']
    #tibia length




#clinical chemistry
Ap4e1_IMPC_CBC_025_001 = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_CBC_025_001']
    #LDL cholesterol

##iron
Ap4e1_Fe = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_CBC_011_001']
##glucose
Ap4e1_Glc = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_CBC_018_001']
##triglycerides
Ap4e1_Tgc = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_CBC_017_001']
##potassium
Ap4e1_K = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_CBC_002_001']
##sodium
Ap4e1_Na = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_CBC_001_001']
##calcium
Ap4e1_Ca = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_CBC_009_001']
##urea
Ap4e1_Urea = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_CBC_011_001']
##total protein
Ap4e1_Prot = Ap4e1[Ap4e1.parameter_stable_id == 'IMPC_CBC_006_001']

#%%
#Nxn

Nxn = pd.read_csv("Nxn.csv")

#body weight (curve)
Nxn_IMPC_BWT_008_001 = Nxn[Nxn.parameter_stable_id == 'IMPC_BWT_008_001']
    #Plot weight curve for gene/sex and each center and compare data

#fertility
Nxn_IMPC_FER_019_001 = Nxn[Nxn.parameter_stable_id == 'IMPC_FER_019_001']
Nxn_IMPC_FER_001_001 = Nxn[Nxn.parameter_stable_id == 'IMPC_FER_001_001']

#viability
Nxn_IMPC_VIA_002  = Nxn[Nxn.parameter_stable_id == 'IMPC_VIA_002']
    #needs to be calculated with formula (excel?)
Nxn_IMPC_VIA_001  = Nxn[Nxn.parameter_stable_id == 'IMPC_VIA_001']

#ECG
Nxn_IMPC_ECG_028_001 = Nxn[Nxn.parameter_stable_id == 'IMPC_ECG_028_001']
    #QT
Nxn_IMPC_ECG_007_001 = Nxn[Nxn.parameter_stable_id == 'IMPC_ECG_007_001']
    #QRS

#IPGTT (Intraperitoneal glucose tolerance test)
Nxn_IMPC_IPG_012_001 = Nxn[Nxn.parameter_stable_id == 'IMPC_IPG_012_001']
    #Area under glucose response curve (AUC)

#X-ray
Nxn_IMPC_XRY_001_001 = Nxn[Nxn.parameter_stable_id == 'IMPC_XRY_001_001']
    #skull shape

#Eye
Nxn_IMPC_EYE_007_001 = Nxn[Nxn.parameter_stable_id == 'IMPC_EYE_007_001']
    #cornea
Nxn_IMPC_EYE_016_001 = Nxn[Nxn.parameter_stable_id == 'IMPC_EYE_0016_001']
    #lens
Nxn_IMPC_EYE_017_001 = Nxn[Nxn.parameter_stable_id == 'IMPC_EYE_017_001']
    #lens opacity

#heart weight
Nxn_IMPC_HWT_008_001 = Nxn[Nxn.parameter_stable_id == 'IMPC_HWT_008_001']
    #heart weight
Nxn_IMPC_HWT_002_001 = Nxn[Nxn.parameter_stable_id == 'IMPC_HWT_002_001']
    #tibia length




#clinical chemistry
Nxn_IMPC_CBC_025_001 = Nxn[Nxn.parameter_stable_id == 'IMPC_CBC_025_001']
    #LDL cholesterol

##iron
Nxn_Fe = Nxn[Nxn.parameter_stable_id == 'IMPC_CBC_011_001']
##glucose
Nxn_Glc = Nxn[Nxn.parameter_stable_id == 'IMPC_CBC_018_001']
##triglycerides
Nxn_Tgc = Nxn[Nxn.parameter_stable_id == 'IMPC_CBC_017_001']
##potassium
Nxn_K = Nxn[Nxn.parameter_stable_id == 'IMPC_CBC_002_001']
##sodium
Nxn_Na = Nxn[Nxn.parameter_stable_id == 'IMPC_CBC_001_001']
##calcium
Nxn_Ca = Nxn[Nxn.parameter_stable_id == 'IMPC_CBC_009_001']
##urea
Nxn_Urea = Nxn[Nxn.parameter_stable_id == 'IMPC_CBC_011_001']
##total protein
Nxn_Prot = Nxn[Nxn.parameter_stable_id == 'IMPC_CBC_006_001']