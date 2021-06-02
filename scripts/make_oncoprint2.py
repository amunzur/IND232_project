#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 16:43:00 2021

@author: amunzur
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np

# load utilities functions to make the oncoprint
exec(open("/groups/wyattgrp/users/amunzur/ind232/scripts/utilities_make_oncoprint.py").read())

########################
# DEFINE VARIABLES
########################

# PATH_sample_ids = "/groups/wyattgrp/users/amunzur/onboarding/data/m1rp_patient_sample_IDs.tsv"
PATH_cn = "/groups/wyattgrp/users/amunzur/ind232/data/ind232_CN.csv"
PATH_muts = "/groups/wyattgrp/users/amunzur/ind232/data/ind232_muts.csv"
PATH_figure = "/groups/wyattgrp/users/amunzur/ind232/figures/ind232_oncoprint.pdf"

patient_id = "" # leave blank (an empty string) to show all patients
filter_ctDNA = False # True means only keep the sample with higher ctDNA content from a patient. Or write False to keep all samples
sort_sample_type = True # gather together EOT and baseline samples. Otherwise no ordering

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
    df_cn = df_cn.drop(columns=['DNA repair defect']) # drop unneeded cols for the oncoprint
    df_cn = pd.melt(df_cn, id_vars = ["Sample", "ctDNA fraction", "Mutation count", "Responder_status"], var_name='Gene', value_name='Copy number')
    
# add two cols for patient id and sample
df_cn[["Patient ID", "Sample type"]] = df_cn["Sample"].str.split(pat = "_", expand = True).rename(columns = {0:"Patient ID", 1:"Sample type"}) # separate sample name from the patient id 
df_muts[["Patient ID", "Sample type"]] = df_muts["Sample"].str.split(pat = "_", expand = True).rename(columns = {0:"Patient ID", 1:"Sample type"}) # separate sample name from the patient id

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
    
    df_cn = df_cn.reset_index(drop = True)
    df_muts = df_muts.reset_index(drop = True)
    
# sort df_cn once more based on ctDNA fraction, from high to low
df_cn = df_cn.sort_values(by = ["ctDNA fraction"])
    
########################
# COLORS and SHAPES
########################
sample_type_dict = {"EOT": "#CC6677", "Baseline": "#DDCC77"}
df_cn["sample_type_color"] = df_cn["Sample type"].map(sample_type_dict)

cn_dict = {-2:'#3f60ac', -1:'#9cc5e9', 0:'#e6e7e8', 1:'#f59496', 2:'#ee2d24'}
df_cn['Color'] = df_cn['Copy number'].map(cn_dict)

mut_dict = {'Missense mutation':'#79B443', 'Splice site mutation':'#FFC907', 'Splice site mutation ': '#FFC907',
               'Stopgain mutation':'#FFC907', 'Frameshift mutation':'#FFC907', 'Frameshift indel':'#5c32a8', 'Non-frameshift indel':'#BD4398',
               'Germline frameshift mutation':'#FFC907', 'Germline stopgain mutation':'#FFC907', 'Germline missense mutation':'#79B443',
               'Intragenic arrangement': "#a3dcf7", "Multiple somatic mutations": "black", "Structural rearrangement": "#FF5733"}

df_muts['Color'] = df_muts['Effect'].map(mut_dict)

shape_dict = {'Missense mutation':'s', 'Splice site mutation':'s', 'Splice site mutation ': 's',
               'Stopgain mutation':'s', 'Frameshift mutation':'s', 'Frameshift indel':'s', 'Non-frameshift indel':'s',
               'Germline frameshift mutation':'*', 'Germline stopgain mutation':'*', 'Germline missense mutation':'*',
               'Intragenic arrangement': "s", "Multiple somatic mutations": "^", "Structural rearrangement": "s"}
df_muts["shapes"] = df_muts["Effect"].map(shape_dict)

########################
# DIVIDE THE DATAFRAMES
########################
[df_cn1, df_cn2] = filter_df_by_col(df_cn, "Responder_status") # not responsive
[df_muts1, df_muts2] = filter_df_by_col(df_muts, "Responder_status") # responsive

[df_cn1_repair, df_cn1_other] = filter_by_genes(df_cn1)
[df_cn2_repair, df_cn2_other] = filter_by_genes(df_cn2)
[df_muts1_repair, df_muts1_other] = filter_by_genes(df_muts1)
[df_muts2_repair, df_muts2_other] = filter_by_genes(df_muts2)

df_counts1_repair = plot_mut_and_cn_counts(df_cn1_repair, df_muts1_repair, drop = True) # not responsive
df_counts1_other = plot_mut_and_cn_counts(df_cn1_other, df_muts1_other, drop = True) # not responsive

df_counts2_repair = plot_mut_and_cn_counts(df_cn2_repair, df_muts2_repair, drop = True) # responsive
df_counts2_other = plot_mut_and_cn_counts(df_cn2_other, df_muts2_other, drop = True) # responsive

# make sure both dfs have the same genes, replace with 0 if exists in one df only
# genes_to_add = df_counts1[~df_counts1.index.isin(df_counts2.index)].dropna().index.to_list()
# df_zero = pd.DataFrame(0, index = genes_to_add, columns = df_counts1.columns) # df of 0s

# df_counts2 = pd.concat([df_counts2, df_zero]) # update by adding a df of zeros to replace missing genes

#Dict to map samples to columns on oncoprint  
samples1 = df_cn1["Sample"].unique().tolist()
samples2 = df_cn2["Sample"].unique().tolist()

repair_genes = ['MSH2', 'MSH6', 'BRCA2', 'CDK12', 'ATM']
other_genes = ['AR', 'SPOP', 'FOXA1', 'TP53', 'RB1', 'PTEN', 'PI3K', 'APC', 'CTNNB1']
gene_pos_repair = {repair_genes[i]: list(range(0,len(repair_genes)))[i] for i in range(len(repair_genes))}
gene_pos_other = {other_genes[i]: list(range(0,len(other_genes)))[i] for i in range(len(other_genes))}

sample_pos1 = {samples1[i]: list(range(0,len(samples1)))[i] for i in range(len(samples1))}
sample_pos2 = {samples2[i]: list(range(0,len(samples2)))[i] for i in range(len(samples2))}

ordered_patients_list1 = [sample.split(sep = "_")[0] for sample in samples1]
ordered_patients_list2 = [sample.split(sep = "_")[0] for sample in samples2]

# order counts dfs based on the gene positions set above
df_counts1_repair = df_counts1_repair.reindex(list(gene_pos_repair.keys()))
df_counts1_other = df_counts1_other.reindex(list(gene_pos_other.keys()))

df_counts2_repair = df_counts2_repair.reindex(list(gene_pos_repair.keys()))
df_counts2_other = df_counts2_other.reindex(list(gene_pos_other.keys()))

# calculate the PERCENTAGE instead of absolute counts
df_counts1_repair = convert_counts_to_percentage(25, 2, 2, df_cn1_repair, df_counts1_repair)
df_counts1_other = convert_counts_to_percentage(25, 2, 2, df_cn1_other, df_counts1_other)

df_counts2_repair = convert_counts_to_percentage(25, 2, 2, df_cn2_repair, df_counts2_repair)
df_counts2_other = convert_counts_to_percentage(25, 2, 2, df_cn2_other, df_counts2_other)


########################
# PREPARE TO PLOT
########################    
fig_width = 8
fig_height = 12

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams["font.size"] = "6"

bar_height = 0.8
bar_width = 0.7

# offset = -(bar_height/3)
offset = -0.4

fig = plt.figure(figsize=(fig_width, fig_height))
gs  = gridspec.GridSpec(nrows=5, ncols=4, height_ratios = [1.8, 1.8, 4, 7, 7], width_ratios = [7, 0.5, 1.3, 0.5], wspace = 0.02, hspace=0.02)
# gs.update(wspace=0.015, hspace=0.05)# set the spacing between axes. 

top1 = fig.add_subplot(gs[0,0]) # ctDNA
mutcount1 = fig.add_subplot(gs[1,0], sharex = top1) # mut count
bottom1_repair = fig.add_subplot(gs[2,0], sharex = mutcount1) # heatmap - repair genes
bottom1_other = fig.add_subplot(gs[3,0], sharex = mutcount1) # heatmap - other genes
count1_repair = fig.add_subplot(gs[2,1], sharey = bottom1_repair)
count1_other = fig.add_subplot(gs[3,1], sharey = bottom1_other, sharex = count1_repair)


top2 = fig.add_subplot(gs[0,2], sharey = top1) # ctDNA
mutcount2 = fig.add_subplot(gs[1,2], sharex = top2, sharey = mutcount1) # mutcount
bottom2_repair = fig.add_subplot(gs[2,2], sharex = mutcount2) #heatmap - repair genes
bottom2_other = fig.add_subplot(gs[3,2], sharex = mutcount2) # heatmap - other genes
count2_repair = fig.add_subplot(gs[2, 3], sharey = bottom2_repair)
count2_other = fig.add_subplot(gs[3, 3], sharey = bottom2_other, sharex = count2_repair)

sub_legend = fig.add_subplot(gs[4,:])

########################
# PLOTTING
########################
# Plot ctDNA fraction in the top subplot
top1.bar(df_cn1["Sample"], df_cn1["ctDNA fraction"], color = "#202020", zorder = 3)
top2.bar(df_cn2["Sample"], df_cn2["ctDNA fraction"], color = "#202020", zorder = 3)

# plot the mutation count 
mutcount1.bar(df_cn1["Sample"], df_cn1["Mutation count"], color = "#202020", zorder = 3)
mutcount2.bar(df_cn2["Sample"], df_cn2["Mutation count"], color = "#202020", zorder = 3)

# plot copy numbers
plot_cn(samples1, repair_genes, df_cn1_repair, bottom1_repair, offset, bar_height, bar_width)
plot_cn(samples1, other_genes, df_cn1_other, bottom1_other, offset, bar_height, bar_width)

plot_cn(samples2, repair_genes, df_cn2_repair, bottom2_repair, offset, bar_height, bar_width)
plot_cn(samples2, other_genes, df_cn2_other, bottom2_other, offset, bar_height, bar_width)

# plot mutations     
plot_muts(sample_pos1, gene_pos_repair, df_muts1_repair, bottom1_repair)
plot_muts(sample_pos1, gene_pos_other, df_muts1_other, bottom1_other)

plot_muts(sample_pos2, gene_pos_repair, df_muts2_repair, bottom2_repair)
plot_muts(sample_pos2, gene_pos_other, df_muts2_other, bottom2_other)

# plot CN and mut count bar graphs 
width = 0.315

cn_repair = dict(zip(gene_pos_repair.keys(), [x - width/2.7 for x in gene_pos_repair.values()]))
mut_repair = dict(zip(gene_pos_repair.keys(), [x - width/1.92 - 0.07 for x in gene_pos_repair.values()]))

cn_other = dict(zip(gene_pos_other.keys(), [x + width/2.7 for x in gene_pos_other.values()]))
mut_other = dict(zip(gene_pos_other.keys(), [x - width/1.92 - 0.07 for x in gene_pos_other.values()]))

count1_repair.barh(list(cn_repair.values()), df_counts1_repair.CN_changes_perc, width, color = "#B0B0B0")
count1_repair.barh(list(mut_repair.values()), df_counts1_repair.Mutation_events_perc, width, color = "#404040")

count1_other.barh(list(cn_other.values()), df_counts1_other.CN_changes_perc, width, color = "#B0B0B0")
count1_other.barh(list(mut_other.values()), df_counts1_other.Mutation_events_perc, width, color = "#404040")

count2_repair.barh(list(cn_repair.values()), df_counts2_repair.CN_changes_perc, width, color = "#B0B0B0")
count2_repair.barh(list(mut_repair.values()), df_counts2_repair.Mutation_events_perc, width, color = "#404040")

count2_other.barh(list(cn_other.values()), df_counts2_other.CN_changes_perc, width, color="#B0B0B0")
count2_other.barh(list(mut_other.values()), df_counts2_other.Mutation_events_perc, width, color="#404040")

########################
# STYLING
########################
for ax in [top2, mutcount2]:
    ax.xaxis.set_visible(False)
    
    # counts bar graphs
for ax in [count1_repair, count2_repair]: 
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(False)
    ax.tick_params(axis='x', length=0)
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.tick_params(axis='y', length=0)
    plt.setp(ax.get_yticklabels(), visible=False)

    # ax.set_xticks([])
    # ax.set_yticks([])

for ax in [count1_other, count2_other]: 
    ax.tick_params(labelrotation = 90, direction = "out", pad = 0.2)
    ax.set_xticks([0, 50, 100])
    ax.set_xticklabels(["0", "50", "100"])
    ax.xaxis.set_tick_params(width = 0.5)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.tick_params(axis='y', length=0)
    plt.setp(ax.get_yticklabels(), visible=False)

top1.set_xticks([])
top1.set_yticks([0, 25, 50, 75, 100])
top1.set_yticklabels(["0", "25", "50", "75", "100"])
top1.set_ylabel("ctDNA %", labelpad=17, rotation = 0, va = 'center')
top1.grid(zorder = 0, linestyle='--', linewidth = 0.5, axis = "y")
top1.tick_params(labelbottom=False)

top2.set_yticks([0, 25, 50, 75, 100])
top2.grid(zorder = 0, linestyle='--', linewidth = 0.5)
plt.setp(top2.get_yticklabels(), visible=False) # remove tick labels from eot plot only

mutcount1.set_xticks([])
mutcount1.set_yticks([0, 25, 50])
mutcount1.set_yticklabels(["0", "25", "50"])
mutcount1.set_ylabel("Total \nmutation \ncount", labelpad=20, rotation = 0, va = 'center')
mutcount1.grid(zorder = 0, linestyle='--', linewidth = 0.5, axis = "y")
mutcount1.tick_params(labelbottom=False)

mutcount2.grid(zorder = 0, linestyle='--', linewidth = 0.5)
plt.setp(mutcount2.get_yticklabels(), visible=False)

bottom1_repair.tick_params(labelrotation = 90, direction = "out", pad = 3)
bottom1_repair.set_yticks(list(range(0, len(repair_genes))))
bottom1_repair.set_yticklabels(repair_genes, rotation = 0, ha = 'right')
bottom1_repair.set_xlim([-1, len(samples1) - 0.5])

bottom2_repair.yaxis.set_visible(False)
bottom2_repair.tick_params(labelrotation = 90, direction = "out", pad = 7)
bottom2_repair.set_xlim([-1, len(samples2) - 0.55])

bottom1_other.tick_params(labelrotation = 90, direction = "out", pad = 3)
bottom1_other.set_yticks(list(range(0, len(other_genes))))
bottom1_other.set_yticklabels(other_genes, rotation = 0, ha = 'right')
bottom1_other.set_xticks(list(range(0, len(ordered_patients_list1))))
bottom1_other.set_xlim([-1, len(samples1) - 0.5])
bottom1_other.set_xticklabels(ordered_patients_list1, ha = "center")

bottom2_other.yaxis.set_visible(False)
bottom2_other.tick_params(labelrotation = 90, direction = "out", pad = 7)
bottom2_other.set_xticks(list(range(0, len(ordered_patients_list2))))
bottom2_other.set_xticklabels(ordered_patients_list2)
bottom2_other.set_xlim([-1, len(samples2) - 0.55])

sub_legend.xaxis.set_visible(False)
sub_legend.yaxis.set_visible(False)
sub_legend.set_facecolor("none")

for ax in [top1, top2, mutcount1, mutcount2, bottom1_repair, bottom2_repair, bottom1_other, bottom2_other, sub_legend]: 
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
   
########################
# LINES
########################

y_position = 0.72
plt.axhline(y = y_position, xmin = 0.01, xmax = bottom1_repair.get_position().x1 + 0.03, color='black', linestyle='-')
plt.axhline(y = y_position, xmin = 0.795, xmax = 0.94, color='black', linestyle='-')

plt.text(0.28, y_position - 0.06, "Best response PD", fontsize = 9)
plt.text(0.84, y_position - 0.06, "PR/SD", fontsize = 9)
    
########################
# LEGEND
########################
    
legend_cn_dict = {"Deep deletion":'#3f60ac', "Deletion":'#9cc5e9', "Neutral":'#e6e7e8', "Gain":'#f59496', "Amplification":'#ee2d24'}

mut_dict = {'Missense':'#79B443', 'Splice site, stopgain, frameshift':'#FFC907', 'Frameshift indel':'#5c32a8', 'Non-frameshift indel':'#BD4398', "Structural rearrangement": "#FF5733"}

mut_dict_shape = {'Somatic':'s', 'Germline':'*', '>2 mutations': '^'}
mut_dict_shape_color = {'Somatic':'#B0B0B0', 'Germline':'#B0B0B0', '>2 mutations': 'black'}


# legend1
handles_cn = []
for key in legend_cn_dict:
    handle = mpatches.Patch(color = legend_cn_dict.get(key), label = key)
    handles_cn.append(handle)

# legend2
handles_muts = []
for key in mut_dict:   
    handle = mpatches.Patch(color = mut_dict.get(key), label = key)
    handles_muts.append(handle)

# legend3
handle_cn = mpatches.Patch(color = "#B0B0B0", label = "Copy number variants")
handle_mut = mpatches.Patch(color = "#404040", label = "Mutations")
handles_counts = [handle_cn, handle_mut]

# legend4
handles_mut_shapes = []
for key in mut_dict_shape:    
    handle = Line2D([0], [0], linestyle = "none", marker = mut_dict_shape.get(key), label = key, markerfacecolor = mut_dict.get(key), color = mut_dict_shape_color.get(key), markersize=5)
    handles_mut_shapes.append(handle)

legend1 = fig.legend(handles=handles_cn, bbox_to_anchor=(0.29, y_position - 0.45), frameon=False, title = "Copy number variants", title_fontsize = 7)
legend2 = fig.legend(handles=handles_muts, bbox_to_anchor=(0.52, y_position - 0.45), frameon=False, title = "Mutations", title_fontsize = 7)
legend3 = fig.legend(handles=handles_counts, bbox_to_anchor=(0.71, y_position - 0.45), frameon=False, title = "Copy number change and \nmutation percentages", title_fontsize = 7)
legend4 = fig.legend(handles=handles_mut_shapes, bbox_to_anchor=(0.42, y_position - 0.53), frameon=False)

# align the legend titlesshapes_dict = {}
legend1._legend_box.align = "left"
legend2.get_title().set_position((-40, 0))
legend3._legend_box.align = "left"


fig.savefig("/groups/wyattgrp/users/amunzur/ind232/figures/fig.pdf", bbox_extra_artists=(), bbox_inches='tight')
plt.show()








