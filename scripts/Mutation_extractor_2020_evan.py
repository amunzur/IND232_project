# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 09:58:49 2018

@author: ewarner

The first part of this script takes a betastasis table (tab separated format) as input
and returns all mutations marked with a * in list format.

The second part is a simple bit that screens a given list of mutations against a list of 
genes of interest - I wrote it to scrape our data for all DDR mutations.

The third part of the script allows for removal of patient mutations listed more than once,
a problem when analyzing mutliple cohorts.
"""
import pandas as pd

PATH_input = "/groups/wyattgrp/users/amunzur/ind232/data/betastasis_tables/rare_germline_variants.tsv"
PATH_output = "/groups/wyattgrp/users/amunzur/ind232/data/betastasis_tables/rare_germline_variants_betastasis_tabulated.csv"

# PATH_input = "/groups/wyattgrp/users/amunzur/ind232/data/betastasis_tables/somatic_mutations.tsv"
# PATH_output = "/groups/wyattgrp/users/amunzur/ind232/data/betastasis_tables/somatic_mutations_betastasis_tabulated.csv"

#ddr list here----------------
ddr_list = ['ATM', 'ATR', 'BRCA1', 'BRCA2', 'CDK12', 'ERCC1', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5',
            'FANCA', 'FANCC', 'FAND2', 'FANCE', 'FANCF', 'FANCG', 'MLH1', 'MSH2', 'MSH6', 'RAD51B',
            'RAD51C']

df= pd.read_csv(PATH_input, sep='\t', header=0)

#Drop columns with uneccessary information
df.drop(['CHROM', 'POSITION', 'REF', 'ALT', 'NOTES'], axis=1, inplace =True)

#remove all commas, they complicate column assignment
df['EFFECT'] = df['EFFECT'].str.replace(',', '')

#reorganize table
df2 = pd.melt(df, id_vars=['GENE', 'EFFECT'], var_name='PATIENT', value_name='FREQUENCY')

#remove entries that don't contain *
df3 =df2[df2['FREQUENCY'].str.contains('*', regex=False)]

####################df3.reset_index(inplace=True)

#Remove all * from output
df3['FREQUENCY'] =df3['FREQUENCY'].str.replace('*', '')

####################df4 = df3.drop(['index'], axis=1)

#Reorder columns
df4 = df3[['PATIENT', 'GENE', 'EFFECT', 'FREQUENCY']]

# Reset index
df4 = df4.reset_index(drop = True)

#Convert to CSV file
df4.to_csv(PATH_output, index=False)

#print ('done part 1')
'''
######################################
# Screen for DDR genes
#####################################

new_df = pd.read_csv('MasterList_Germline.csv', sep=',', header=0)

new_df2 = new_df[new_df['GENE'].isin(ddr_list)]

new_df2.to_csv('list_DDR_germline.csv', sep ='\t', index=False)

print ('done part 2')


#####################################
# Searching for overlapping sample ID's
#####################################

data_f = pd.read_csv("MasterList_Germline_wGU_m1.tsv", sep='\t', header=0)


data_f['INFO'] = data_f['GENE'] + data_f['EFFECT']

data_f.drop(['GENE', 'EFFECT', 'FREQUENCY'], axis=1, inplace=True)

mask=data_f['INFO'].duplicated(keep=False)

df_mask = data_f[mask]

df_mask.to_csv('wtf.csv', sep='\t', index=False)
print(mask)
'''





