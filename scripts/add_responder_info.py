#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 21 10:56:11 2021

@author: amunzur
"""

import pandas as pd 

PATH_resp = "/groups/wyattgrp/users/amunzur/ind232/data/2020-09-03_I232_Responders_Jean.xlsx"
PATH_muts = "/groups/wyattgrp/users/amunzur/ind232/data/ind232_muts.csv"
PATH_cn = "/groups/wyattgrp/users/amunzur/ind232/data/ind232_CN.csv"

df_resp = pd.read_excel(PATH_resp)
df_muts = pd.read_csv(PATH_muts)
df_cn = pd.read_csv(PATH_cn)

# add a patients column
df_cn = pd.concat([df_cn["Sample"].str.split("_", expand = True), df_cn.reset_index(drop = True)], axis = 1) # choose muts repeated at least twice, split on "and" and expand into separate cols
df_cn = df_cn.rename(columns = {0: "Patient ID", 1: "Type"})

df_muts = pd.concat([df_muts["Sample"].str.split("_", expand = True), df_muts.reset_index(drop = True)], axis = 1) # choose muts repeated at least twice, split on "and" and expand into separate cols
df_muts = df_muts.rename(columns = {0: "Patient ID", 1: "Type"})

# make sure the syntax is the same in between dataframes
responders = df_resp["Patient"].unique()
responders = [x[:4] + "-" + x[4:] for x in responders] # 

# add a new col for responders 
df_cn["Responder_status"] = df_cn["Patient ID"].isin(responders)
df_muts["Responder_status"] = df_muts["Patient ID"].isin(responders)

# replace true with responder and false with non responder 
df_cn["Responder_status"] = df_cn["Responder_status"].replace({True: 'Responsive', False: 'Not responsive'})
df_muts["Responder_status"] = df_muts["Responder_status"].replace({True: 'Responsive', False: 'Not responsive'})

# order 
df_cn = df_cn.sort_values(by = ["Responder_status"])
df_muts = df_muts.sort_values(by = ["Responder_status"])

# drop the new cols we made 
df_cn = df_cn.drop(["Patient ID", "Type"], axis = 1)
df_muts = df_muts.drop(["Patient ID", "Type"], axis = 1)

# move the responder col to the first 
cols = list(df_cn.columns)
cols = [cols[-1]] + cols[:-1]
df_cn = df_cn[cols]

cols = list(df_muts.columns)
cols = [cols[-1]] + cols[:-1]
df_muts = df_muts[cols]

# save the updated dataframes
df_cn.to_csv(PATH_cn, index = False)
df_muts.to_csv(PATH_muts, index = False)