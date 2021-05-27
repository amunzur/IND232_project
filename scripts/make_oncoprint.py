#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 16:43:00 2021

@author: amunzur
"""

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap


########################
# DEFINE VARIABLES
########################

# PATH_sample_ids = "/groups/wyattgrp/users/amunzur/onboarding/data/m1rp_patient_sample_IDs.tsv"
PATH_cn = "/groups/wyattgrp/users/amunzur/ind232/data/ind232_CN.csv"
PATH_muts = "/groups/wyattgrp/users/amunzur/ind232/data/ind232_muts.csv"
PATH_figure = "/groups/wyattgrp/users/amunzur/ind232/figures/ind232_oncoprint.pdf"

patient_id = "" # leave blank (an empty string) to show all patients
filter_ctDNA = True # True means only keep the sample with higher ctDNA content from a patient. Or write False to keep all samples
sort_sample_type = True # gather together EOT and baseline samples. Otherwise no ordering

fig_width = 5
fig_height = 4

########################
# DEFINE FUNCTIONS
########################

def filter_df_ctDNA(df):
    '''
    Filter the data frame based on ctDNA content. From baseline or EOT samples only keep the one
    with higher ctDNA content.
    '''
    # from the same patient only keep the higher of the two samples for both dfs 
    idx = list(df.groupby(['Patient ID'], sort=False)['ctDNA fraction'].transform(max) == df["ctDNA fraction"])
    df = df[idx].reset_index(drop = True)
    
    return(df)

########################
# READ FILES
########################

# read the dataframes based on type
if PATH_cn.split(".")[1] == "tsv": 
    df_cn = pd.read_csv(PATH_cn, sep = "\t")
    df_muts = pd.read_csv(PATH_muts, sep = "\t")

elif PATH_cn.split(".")[1] == "csv":
    df_cn = pd.read_csv(PATH_cn)
    df_muts = pd.read_csv(PATH_muts)

# if patient id is given, filter for that 
if patient_id:
    df_cn = df_cn.loc[df_cn['Patient ID'] == patient_id]
    df_muts = df_muts.loc[df_muts['Patient ID'] == patient_id]
else:
    pass

########################
# MODIFY FILES
########################

# make sure the genes are in the rows and samples are columns 
if df_cn.index.name is None: # if the index hasn't been set
    df_cn = df_cn.drop(columns=['Unnamed: 0', 'DNA repair defect']) # drop unneeded cols for the oncoprint
    df_cn = pd.melt(df_cn, id_vars = ["Sample", "ctDNA fraction", "Mutation count"], var_name='Gene', value_name='Copy number')
    
for df in [df_cn, df_muts]: 
    df[["Patient ID", "Sample type"]] = df["Sample"].str.split(pat = "_", expand = True).rename(columns = {0:"Patient ID", 1:"Sample type"}) # separate sample name from the patient id 

df_cn["ctDNA fraction"] = df_cn["ctDNA fraction"] * 100 # convert fraction to percentage 

# make a new df for sample_type colors 
df_type = df_cn.drop_duplicates(subset = ["Sample"])

# filtering
if filter_ctDNA == True:
    df_cn = filter_df_ctDNA(df_cn) # filter based on ctDNA content
    df_muts = df_muts[df_muts["Sample"].isin(df_cn["Sample"])] # filter df_muts to keep samples from df_cn
else: 
    pass

if sort_sample_type == True:
    
    # ordering based on sample type 
    df_cn = df_cn.sort_values(by = "Sample type")
    df_muts = df_muts.sort_values(by = "Sample type")
    
    # order based on ctDNA content, within the groups
    df_cn = df_cn.groupby("Sample type").apply(pd.DataFrame.sort_values, 'ctDNA fraction')
    # df_muts = df_muts.groupby("Sample type").apply(pd.DataFrame.sort_values, 'ctDNA fraction')
    
########################
# MAPPING
########################

# COLORS
sample_type_dict = {"EOT": "#CC6677", "Baseline": "#DDCC77"}
df_cn["sample_type_color"] = df_cn["Sample type"].map(sample_type_dict)

cn_dict = {-2:'#3f60ac', -1:'#9cc5e9', 0:'#e6e7e8', 1:'#f59496', 2:'#ee2d24'}
df_cn['Color'] = df_cn['Copy number'].map(cn_dict)

mut_dict = {'Missense mutation':'#79B443', 'Frameshift mutation':'#BD4398', 'Splice site mutation':'#FFC907',
               'Stopgain mutation':'#FFC907', 'Germline frameshift mutation':'#8c69ff', 'Germline stopgain mutation':'#FF5733'}
df_muts['Color'] = df_muts['Effect'].map(mut_dict)

########################
# PREPARE TO PLOT
########################
bar_height = 0.7
bar_width = 0.7

samples = df_cn["Sample"].unique().tolist()
genes = df_cn['Gene'].unique().tolist()

ordered_patients_list = [sample.split(sep = "_")[0] for sample in samples]

# offset = -(bar_height/3)
offset = -0.4

fig = plt.figure(figsize=(fig_width, fig_height))
gs  = gridspec.GridSpec(nrows=4, ncols=1, height_ratios = [2, 2, 0.4, 15], hspace = 0.1, wspace = 0)
# gs.update(wspace=0.015, hspace=0.05)# set the spacing between axes. 

sub_top = fig.add_subplot(gs[0,0]) # ctDNA
sub_mutcount = fig.add_subplot(gs[1,0]) # mut count
sub_type = fig.add_subplot(gs[2,0]) # sample type 
sub_bottom = fig.add_subplot(gs[3,0]) # heatmap

#Dict to map genes to row on oncoprint               
gene_pos = {genes[i]: list(range(0,len(genes)))[i] for i in range(len(genes))}
sample_pos = {samples[i]: list(range(0,len(samples)))[i] for i in range(len(samples))}

########################
# PLOTTING
########################
# Plot ctDNA fraction in the top subplot
sub_top.bar(df_cn["Sample"], df_cn["ctDNA fraction"], color = "#202020", zorder = 3)

# plot the sample types (EOT or baseline)
sub_type.bar(x = df_cn["Sample"], height = bar_height, color = df_cn["sample_type_color"], width=0.8, edgecolor=None, linewidth=0, zorder = 3)

# plot the mutation count 
sub_mutcount.bar(df_cn["Sample"], df_cn["Mutation count"], color = "#202020", zorder = 3)

# Plot copy numbers
for sample in samples:
    bottom = offset 
    for gene in genes:
        row = df_cn.loc[(df_cn['Gene'] == gene) & (df_cn['Sample'] == sample)]
        color = row['Color'].values[0]

        sub_bottom.bar(sample, bar_height, bottom = bottom, color = color, zorder = 10, width = bar_width * 1.2)

        bottom += 1

# Plot mutations
for i, row in df_muts.iterrows():
    sample = row['Sample']
    mut_type = row['Effect']
    gene = row['Gene']
    color = row['Color']
    marker_type = "s"
    
    if "germ" in mut_type: marker_type = "*"    
    
    # check if there is another mutation in the same sample/gene combination
    x = [sample, gene] == df_muts[["Sample", "Gene"]]
    if sum(x.all(axis = 1)) == 1: # ONE mutation only in the same sample and gene combination        
        sub_bottom.scatter(x = sample_pos[sample], y = gene_pos[gene], c = color, s = 3, marker = marker_type, zorder = 100)

    elif sum(x.all(axis = 1)) == 2: # TWO mutations in the same sample and gene combination
        sub_bottom.scatter(x = sample_pos[sample] + 0.08, y = gene_pos[gene] + 0.08, c = color, s = 3, marker = marker_type, zorder = 100)
        sub_bottom.scatter(x = sample_pos[sample] - 0.04, y = gene_pos[gene] - 0.08, c = color, s = 3, marker = marker_type, zorder = 100)
    
    elif sum(x.all(axis = 1)) == 3: # THREE mutations 
        sub_bottom.scatter(x = sample_pos[sample], y = gene_pos[gene], c = "black", s = 6, marker = "^", zorder = 100)
        
########################
# STYLING
########################
sub_top.set_xticks([])
sub_top.set_yticks([0, 50, 100])
sub_top.set_yticklabels(["0", "50", "100"])
sub_top.set_ylabel("ctDNA %", labelpad=17, rotation = 0, va = 'center')
sub_top.grid(zorder = 0, linestyle='--', linewidth = 0.5)

sub_mutcount.set_xticks([])
sub_mutcount.set_yticks([0, 25, 50])
sub_mutcount.set_yticklabels(["0", "25", "50"])
sub_mutcount.set_ylabel("Mutation \n count", labelpad=20, rotation = 0, va = 'center')
sub_mutcount.grid(zorder = 0, linestyle='--', linewidth = 0.5)

sub_type.xaxis.set_visible(False)
# sub_type.yaxis.set_visible(False)
sub_mutcount.set_xticks([])
sub_type.set_yticks([])
sub_type.set_ylabel("Sample type", labelpad=30, rotation = 0, va = 'center')

sub_bottom.tick_params(labelrotation = 90, direction = "out")
sub_bottom.set_yticks(list(range(0, len(genes))))
sub_bottom.set_yticklabels(genes, rotation = 0)
sub_bottom.set_xticks(list(range(0, len(ordered_patients_list))))
sub_bottom.set_xticklabels(ordered_patients_list)

for ax in [sub_top, sub_mutcount, sub_type, sub_bottom]: 
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.tick_params(axis='both', which='major')

    ax.set_xlim([-1, len(samples)])

# set up font size 
font = {'family' : 'normal', 'weight' : 'normal', 'size'   : 6}
plt.rc('font', **font)
# plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.1)
fig.tight_layout(pad=2)
fig.savefig("/groups/wyattgrp/users/amunzur/ind232/figures/fig.pdf", bbox_extra_artists=(), bbox_inches='tight')









