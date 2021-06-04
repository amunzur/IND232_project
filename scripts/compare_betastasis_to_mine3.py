#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 12:11:42 2021

@author: amunzur

The purpose of this script is: 
    1. Reformat the betastasis files in such a way that they are in the same format as mine 
    2. Go through both files to see if we missed any mutations 
"""

import pandas as pd

PATH_beta = "/groups/wyattgrp/users/amunzur/ind232/data/betastasis_tables/somatic_mutations_betastasis_tabulated.csv"
PATH_muts = "/groups/wyattgrp/users/amunzur/ind232/data/ind232_muts.csv"

beta = pd.read_csv(PATH_beta)
muts = pd.read_csv(PATH_muts)

################################
# PART 1. reformat already formatted betastasis
################################
beta.PATIENT = beta.PATIENT.str.replace("cfDNA", "EOT")

for index, row in beta.iterrows():
    
    # deal with patient naming
    x = beta.PATIENT.str.split("-")[index]
    beta.at[index, "PATIENT"] = x[1] + "-" + x[2] + "_" + x[3]
    
    # deal with mut effects naming 
    y = beta.EFFECT.str.split(" ")[index][0]
    beta.at[index, "EFFECT"] = y
    
for index, row in beta.iterrows():
    
    effect = beta.at[index, "EFFECT"]
    if "Non-frameshift" in effect: 
        beta.at[index, "EFFECT"] = effect + " indel"
    elif "Splice" in effect: 
        beta.at[index, "EFFECT"] = effect + " site"
        
beta.PATIENT = beta.PATIENT.str.replace("00", "000")
beta.PATIENT = beta.PATIENT.str.replace("010", "0010")



# del beta["FREQUENCY"]

# only keep the genes we are interested in 
genes_list = list(muts.Gene.unique())
beta = beta[beta["GENE"].isin(genes_list)]


################################
# PART 2. reformat mine, just minor fixes
################################
del muts["Responder_status"]
muts.columns = ["PATIENT", "GENE", "EFFECT"]
muts.EFFECT = muts.EFFECT.str.replace(" mutation", "")
muts.EFFECT = muts.EFFECT.str.rstrip()


################################
# PART 3. compare files
################################
df = pd.merge(beta, muts, how="outer", on = ['PATIENT', 'GENE', 'EFFECT'], indicator=True)
df = df[df['_merge'] == 'left_only']

df = df[df['PATIENT'].isin(muts.PATIENT)]


df = df.merge(muts.PATIENT, how = "inner", on = "PATIENT")

df.to_csv("/groups/wyattgrp/users/amunzur/ind232/data/betastasis_only_mutations.csv")






