#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 21 11:37:12 2021

@author: amunzur
"""


def filter_df_ctDNA(df):
    '''
    Filter the data frame based on ctDNA content. From baseline or EOT samples only keep the one
    with higher ctDNA content.
    '''
    # from the same patient only keep the higher of the two samples for both dfs 
    idx = list(df.groupby(['Patient ID'], sort=False)['ctDNA fraction'].transform(max) == df["ctDNA fraction"])
    df = df[idx].reset_index(drop = True)
    
    return(df)

def filter_df_type(df, key, colname):
    '''
    Filter the data frame based on sample type content.
    '''
    idx = list(df[colname] == key)
    new_df = df[idx]
    
    return(new_df)

def plot_cn(samples, genes, df_cn, ax, offset, bar_height, bar_width):
    '''
    # Plot copy number variations one by one, creating a heatmap.
    '''
    
    for sample in samples:
        bottom = offset 
        for gene in genes:
            row = df_cn.loc[(df_cn['Gene'] == gene) & (df_cn['Sample'] == sample)]
            color = row['Color'].values[0]
            ax.bar(sample, bar_height, bottom = bottom, color = color, zorder = 10, width = bar_width * 1.2)
            bottom += 1
        
def plot_muts(sample_pos, gene_pos, df_muts, ax, offset_how):
    '''
    # Add mutations to the CN heatmap.
    '''
    
    # map a dict
    lw_dict = {"s": 0, "*": 1, "^": 1}
    df_muts['lw_col'] = df_muts['shapes'].map(lw_dict)

    marker_size = 30
    
    for i, row in df_muts.iterrows():
        sample = row['Sample']
        gene = row['Gene']
        color = row["Color"]

        x = [sample, gene] == df_muts[["Sample", "Gene"]]
        df = df_muts[x.all(axis = 1)].reset_index(drop = True) # subset to that gene and sample only 
            
        if df["Effect"].str.contains("Multiple somatic mutations").any():
            if len(df) == 1: 
                ax.scatter(x = sample_pos[sample], y = gene_pos[gene] + 0.02, c = "black", s = marker_size + 3, marker = "^", zorder = 100, lw = 0)
        else:
            if len(df) == 1:
                if "Germline" in row["Effect"]:
                    ax.scatter(x = sample_pos[sample], y = gene_pos[gene] + 0.02, c = color, s = marker_size, marker = str(df["shapes"][0]), zorder = 100, lw = 1)
                else: 
                    ax.scatter(x = sample_pos[sample], y = gene_pos[gene] + 0.02, c = color, s = marker_size, marker = str(df["shapes"][0]), zorder = 100, lw = 0)

            elif len(df) == 2: # two mutations
                if offset_how == "vertical": 
                    ax.scatter(x = sample_pos[sample], y = gene_pos[gene] - 0.17, c = df.loc[0, "Color"], s = marker_size, marker = df.loc[0, "shapes"], zorder = 100, lw = df.loc[0, "lw_col"])
                    ax.scatter(x = sample_pos[sample], y = gene_pos[gene] + 0.17, c = df.loc[1, "Color"], s = marker_size, marker = df.loc[1, "shapes"], zorder = 100, lw = df.loc[1, "lw_col"])
                else: 
                    ax.scatter(x = sample_pos[sample] - 0.07, y = gene_pos[gene] - 0.03, c = df.loc[0, "Color"], s = marker_size, marker = df.loc[0, "shapes"], zorder = 100, lw = df.loc[0, "lw_col"])
                    ax.scatter(x = sample_pos[sample] + 0.07, y = gene_pos[gene] + 0.05, c = df.loc[1, "Color"], s = marker_size, marker = df.loc[1, "shapes"], zorder = 100, lw = df.loc[1, "lw_col"])
                
            elif len(df) >= 2: # three mutations3
                ax.scatter(x = sample_pos[sample], y = gene_pos[gene] + 0.02, c = "black", s = marker_size + 3, marker = "^", zorder = 100, lw = 0)
            
def plot_mut_and_cn_counts(df_cn, df_muts, drop): 
    '''
    Given a df_cn and df_muts data frames, plot the number of CN changes and mutation events in vertical bar graphs. 
    '''
    # get a list of genes from df_cn 
    genes = df_cn["Gene"].unique()
    
    # from df cn, drop rows if CN is 0
    df_cn = df_cn.drop(df_cn[df_cn["Copy number"] == 0].index)
    df_cn = df_cn.groupby(['Gene']).size()
    
    if drop == True:
        df_muts = df_muts.drop_duplicates(['Sample','Gene'],keep= 'last')
    else:
        pass

    df_muts = df_muts.groupby(['Gene']).size()
        
    # combine the dfs into one 
    combined = pd.concat([df_cn, df_muts], axis=1)
    combined.columns = ["CN_changes", "Mutation_events"]
    
    combined = combined.fillna(0)
    
    #add the genes we lost above 
    for gene in genes: 
        if gene not in combined.index:      
            row = pd.DataFrame({'Gene': [gene],
                   'CN_changes': [0],
                   'Mutation_events': [0]}).set_index("Gene")
            
            combined = pd.concat([combined, row])
        
    return combined

def filter_df_by_col(df, colname):
    '''
    Given a df and a col, return separated dfs with the unique names in the given col. 
    '''
    keys = df[colname].unique()
    
    df_list = []
    for key in keys: 
        idx = list(df[colname] == key)
        filtered_df = df[idx]
        df_list.append(filtered_df)
        
    return df_list
    
    
def convert_counts_to_percentage(CN_filter, mut_filter, AR_filter, df_cn, df_counts):
    '''
    Convert mut counts and CN changes to percentages instead of absolute counts.
    '''
    
    dups = df_cn.drop_duplicates(subset = "Sample") # only keep one patient, drop duplicates
    CN_filt = len(dups.loc[dups["ctDNA fraction"] > CN_filter])
    mut_filt = len(dups.loc[dups["ctDNA fraction"] > mut_filter])
        
    # add to the dfs 
    df_counts["CN_changes_perc"] = round(df_counts.CN_changes/CN_filt * 100)
    df_counts["Mutation_events_perc"] = round(df_counts.CN_changes/mut_filt * 100)
        
    # deal with AR, if it's in the df indices
    if df_counts.index.str.contains("AR").any(): 
        AR_filt = len(dups.loc[dups["ctDNA fraction"] > AR_filter])
        df_counts.loc["AR", "CN_changes_perc"] = round((df_counts.loc["AR", "CN_changes"] / AR_filt) * 100)
        df_counts.loc["AR", "Mutation_events_perc"] = round((df_counts.loc["AR", "Mutation_events"] / AR_filt) * 100)
    
    return df_counts
    
    
def filter_by_genes(df, repair, other): 
    '''
    Given a df_cn and muts data frames, separate them based on genes. One df will be DNA repair genes, other will be the rest. 
    Two data frames are returned. 
    '''
    
    idx_repair = df["Gene"].isin(repair)
    idx_other = df["Gene"].isin(other)
    
    df_repair = df[idx_repair]
    df_other = df[idx_other]
        
    return[df_repair, df_other]

    
    

    
    
    
    
    
    
    
    
    
    
    
    