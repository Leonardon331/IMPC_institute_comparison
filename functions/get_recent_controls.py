# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 21:33:30 2023

@author: Leonard Borchardt
"""

import os
import pandas as pd
from urllib.error import HTTPError, URLError
from datetime import datetime, timedelta

#set wd
os.chdir(r'C:\Users\Leonard Borchardt\OneDrive\Studium\Internships\Phenogenomics\Bachelor thesis')
#%%

def get_recent_controls(days: int = 31, months_ago: int = False):
    '''
    gets recent or from certain month back in time IMPC control data for past month via SOLR API and saves it as a csv file

    Parameters
    ----------
    days: int, amount of days for getting control data, default: 31, optional
    months_ago: int, get data from certain month back in time, optional


    Returns
    -------
    None.

    '''
    
    save = r'./controls'
    # check whether directory already exists, otherwise create new folder
    if not os.path.exists(save):
        os.mkdir(save)
    
    
    # Get today's date
    today = datetime.utcnow()
    
    if months_ago:
        # Calculate the date certain month ago
        day_months_ago = today - timedelta(days=months_ago*30)
        
        # Set the time to 00:00:00
        midnight = day_months_ago.replace(hour=0, minute=0, second=0, microsecond=0)
        
    else:
        midnight = today.replace(hour=0, minute=0, second=0, microsecond=0)
    
    # Generate a list of dates for the past 10 k days at 00:00:00
    past_dates_midnight = [midnight - timedelta(days=i) for i in range(3653)]
    
    # Format the dates in ISO 8601 format
    iso_formatted_dates = [date.strftime('%Y-%m-%dT%H:%M:%SZ') for date in past_dates_midnight]
    
    dates = iso_formatted_dates
    #print(dates)
    

    ### get recent control data ###
    # get control data for the most recent timeframe_days (default 31 days) in the past year
    WT = []
    #centers which much data: get time dependent data
    for center in ["TCP", "JAX", "HMGU", "CCP-IMG", "MRC Harwell", "MARC", "KMPC", "RBRC", "ICS", "WTSI", "UC Davis", "BCM", "UCD", "RIKEN BRC"]:
        count = 0
        for date in dates:
                try:
                    #url = f'https://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=biological_sample_group:control%20AND%20date_of_birth:"{date}"%20AND%20phenotyping_center:"{center}"&fl=phenotyping_center,parameter_name,parameter_stable_id,sex,data_point,time_point,biological_sample_group,colony_id,specimen_id,project_name,date_of_birth,age_in_weeks,age_in_days,developmental_stage_name&wt=csv&rows=100000000'
                    url = f'https://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=biological_sample_group:control%20AND%20date_of_birth:"{date}"%20AND%20phenotyping_center:"{center}"&wt=csv&rows=100000000'
                    #centers with space (e.g. UC Davis) needs to be in double quotation marks "" for url, date needs to be surrounded by "", too
                    if " " in url:
                        print("replace space in url")
                        url = url.replace(" ", "%20")
            
                    wt = pd.read_csv(url)
                    
                    if wt.empty:
                        print(f"no control data at day {date} at {center}")
                    
                    else:
                        WT.append(wt)
                        print(f"WT data for {date} at {center} successfully collected!")
                        count+=1
                        if count == days:
                            print(f"Got all control data at {center} for the nearest past {days} control birthdays!")
                            print("-----------------------------------------------------------------------")
                            break
                        
                except (HTTPError, URLError):   #, TimeoutError
                #check for multiple errors with tuple
                    print("url error")
        
    #centers with few data: get time independent data
    for center in ["Monterotondo", "SEAT", "CDTA", "CIPHE", "NARLabs", "Monash", "Oulu", "CAM-SU GRC", "VETMEDUNI", "Monterotondo R&D", "CNB"]: 
        try:
            #url = f'https://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=biological_sample_group:control%20AND%20phenotyping_center:"{center}"&fl=phenotyping_center,parameter_name,parameter_stable_id,sex,data_point,time_point,biological_sample_group,colony_id,specimen_id,project_name,date_of_birth,age_in_weeks,age_in_days,developmental_stage_name&wt=csv&rows=100000000'
            url = f'https://www.ebi.ac.uk/mi/impc/solr/experiment/select?q=biological_sample_group:control%20AND%20phenotyping_center:"{center}"&wt=csv&rows=100000000'
            if " " in url:
                print("replace space in url")
                url = url.replace(" ", "%20")

            wt = pd.read_csv(url)
            
            if wt.empty:
                print(f"no control data at {center}")
            
            else:
                WT.append(wt)
                print(f"All WT data for {center} successfully collected!")
                print("-----------------------------------------------------------------------")
                     
        
        except (HTTPError, URLError):   #, TimeoutError
        #check for multiple errors with tuple
            print("url error")
       
        
    #concat list of data frames (wt_list) to one data frame (WT)
    try:
        WT = pd.concat(WT, ignore_index=True)
        WT = WT.sort_values(by=["procedure_name","parameter_stable_id", "sex","phenotyping_center"]).reset_index(drop=True)
        WT.to_csv(f'{save}/recent_controls.csv', index=False)
        
        #globals()[r"recent_controls"] = WT
        #automatic variable assignment for the specific controls
        
    except ValueError:
        #if WT is still emty (pd.concat will raise ValueError)
        print("no control data at all")
            
   
#%%
get_recent_controls(days=60,months_ago=5)