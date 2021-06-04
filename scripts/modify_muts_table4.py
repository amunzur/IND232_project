#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 15:59:30 2021

@author: amunzur
"""
import pandas as pd 
import numpy as np 

PATH_muts = "/groups/wyattgrp/users/amunzur/ind232/data/ind232_muts.csv"
PATH_new_muts = ""

df_muts = pd.read_csv(PATH_muts)
df_cn = pd.read_csv(PATH_muts)

# remove the SNV in CAVK-0005
df_muts = df_muts.drop([42]) # CAVK-0005 - PTEN
df_muts = df_muts.drop([10]) # CAVA-0003 - FOXA1
df_muts = df_muts.drop([91]) # CAMP-0007 - APC
df_muts = df_muts.drop([44, 45]) # CAVK-0005 - TP53, SPOP

df_muts.at[100, "Color"] = "#AED28E"
df_muts.at[26, "Color"] = "#AED28E"

to_add = pd.DataFrame([["Not responsive", "CAVK-0005_Baseline", "APC", "Stopgain mutation", "CAVK-0005", "Baseline", "#FFC907", "s"], 
                         ["Not responsive", "CAMP-0003_Baseline", "APC", "Missense mutation", "CAMP-0003", "Baseline", "#79B443", "s"], 
                         ["Not responsive", "CAMP-0003_Baseline", "CTNNB1", "Missense mutation", "CAMP-0003", "Baseline", "#79B443", "s"], 
                         ["Not responsive", "CAKO-0005_EOT", "CTNNB1", "Missense mutation", "CAKO-0005", "EOT", "#79B443", "s"]])

to_add.columns = df_muts.columns
df_muts = df_muts.append(to_add, ignore_index = True)    

to_add = pd.DataFrame([["CAVA-0008_Baseline", 59, 8, "Not responsive", "MSH2", -1, "CAVA-0008", "Baseline", "#DDCC77", "#e6e7e8"]])
to_add.columns = df_cn.columns
df_cn = df_cn.append(to_add, ignore_index = True)    
