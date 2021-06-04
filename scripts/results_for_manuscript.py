#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 14:03:26 2021

@author: amunzur
"""
import pandas as pd 
import statistics as stats
exec(open("/groups/wyattgrp/users/amunzur/ind232/scripts/utilities_make_oncoprint.py").read())


patient_id = "" # leave blank (an empty string) to show all patients
filter_ctDNA = True # True means only keep the sample with higher ctDNA content from a patient. Or write False to keep all samples
sort_sample_type = False # gather together EOT and baseline samples.


PATH_cn = "/groups/wyattgrp/users/amunzur/ind232/data/ind232_CN.csv"
df_cn = pd.read_csv(PATH_cn)

if df_cn.index.name is None: # if the index hasn't been set
    df_cn = df_cn.drop(columns=['DNA repair defect']) # drop unneeded cols for the oncoprint
    df_cn = pd.melt(df_cn, id_vars = ["Sample", "ctDNA fraction", "Mutation count", "Responder_status"], var_name='Gene', value_name='Copy number')
    
# add two cols for patient id and sample
df_cn[["Patient ID", "Sample type"]] = df_cn["Sample"].str.split(pat = "_", expand = True).rename(columns = {0:"Patient ID", 1:"Sample type"}) # separate sample name from the patient id 

df_cn["ctDNA fraction"] = df_cn["ctDNA fraction"] * 100 # convert fraction to percentage 

# make a new df for sample_type colors 
df_type = df_cn.drop_duplicates(subset = ["Sample"])

# filtering
if filter_ctDNA == True:
    df_cn = filter_df_ctDNA(df_cn) # filter based on ctDNA content
else: 
    pass

if sort_sample_type == True:
    
    # ordering based on sample type 
    df_cn = df_cn.sort_values(by = "Sample type")
    
    # order based on ctDNA content, within the groups
    df_cn = df_cn.groupby("Sample type").apply(pd.DataFrame.sort_values, 'ctDNA fraction')
    df_cn = df_cn.reset_index(drop = True)
    
# drop a patient 
df_cn = df_cn[df_cn.Sample != "CALM-0002_EOT"]
df = df_cn.drop(['ctDNA fraction', 'Mutation count', 'Responder_status',
       'Copy number', 'Patient ID', 'Sample type'], axis = 1)
df = df.pivot(index = 'Sample', columns='Gene')
idx = (df.index).to_list()

#########################################
# NUMBER OF EOT AND BASELINE SAMPLES USED
#########################################

EOT_number = sum("EOT" in x for x in idx)
Baseline_number = sum("Baseline" in x for x in idx)

#########################################
# MEAN & MEDIAN BASELINE ctDNA FRACTION
#########################################
baseline_df = df_cn[df_cn["Sample type"] != "EOT"] # only keep the baseline samples 
baseline_df = baseline_df.drop_duplicates(subset = ["Sample"])
BL_ctDNA = baseline_df["ctDNA fraction"].to_list()
BL_ctDNA_mean = stats.mean(BL_ctDNA)
BL_ctDNA_median = stats.median(BL_ctDNA)

#########################################
# NUMBER OF EOT AND BASELINE SAMPLES USED
#########################################

print("Number of EOT samples:", EOT_number)
print("Number of baseline samples:", Baseline_number)
print("Total number of samples:", EOT_number + Baseline_number)

print("Baseline samples median ctDNA percentage:", BL_ctDNA_median)
print("Baseline samples mean ctDNA percentage:", BL_ctDNA_mean)















