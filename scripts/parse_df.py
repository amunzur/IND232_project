import numpy as np
import pandas as pd 

#######################
# DEFINE FUNCTIONS
#######################

def cleanup_gene(gene_name, df, colname):
    '''
    Search for the given df and colname for the mutations that begin with the gene name,
    then remove the gene name and retain only the mutation effect.
    '''
    gene_name = gene_name.ljust(len(gene_name) + 1) # add a space to the end of the string
    idx = df_muts[colname].str.startswith(gene_name)
    df.loc[idx, colname] = df.loc[idx, colname].str.split(gene_name).str[1]
    
    return(df)

def remove_locations(to_search, to_replace, df, colname):
    '''
    Given a df and mutations ids, searh for the given mutation type and if it contains the 
    location of the mutation, starting with (p...) remove the location information
    '''
    df_muts.loc[df_muts[colname].astype(str).str.contains(to_search, na = False, case = False), colname] = to_replace
    
    return(df)


#######################
# MODIFY THE MAIN DF
#######################

PATH_to_df = "/groups/wyattgrp/users/amunzur/ind232/data/IND232_cfDNAresults_Aug2019.csv"

df_main = pd.read_csv(PATH_to_df)

# minor modifications to rename and drop unwanted cols 
df_main = df_main.loc[df_main['ctDNA detected'] == "yes"] # remove patients that didn't have any ctDNA detected 
df_main = df_main.sort_values(by = ["Study ID"]) # sort based on patient id 

df_main['Sample'] = df_main[df_main.columns[[0, 1]]].apply(lambda x: '_'.join(x.astype(str)), axis=1) # merge study id and sample cols into one

df_main = df_main.drop(["Study ID", "ctDNA detected", "Notes and other less common alterations"], axis = 1) # drop some unwanted cols
nameslist = ['Sample', 'ctDNA fraction', 'Mutation count', 'DNA repair defect', 'BRCA2', 'ATM', 'CDK12','AR', 'AR_muts', 'PTEN', 'PI3K/AKT', 'RB1','TP53', 'APC_or_CTNNB1', 'SPOP', 'FOXA1']
df_main.columns = nameslist

# clearing up some confusion about CN Loss types
# idx = df_main[df_main["Sample"] == "CAMP-0004_EOT"].index.values
# df_main.loc[idx, "PTEN"] = "CN loss (monoallelic)"

# divide a col into separate cols based on strings 
df_main["APC"] = df_main["APC_or_CTNNB1"]
df_main["CTNNB1"] = df_main["APC_or_CTNNB1"]
del df_main["APC_or_CTNNB1"]

df_main.loc[df_main["APC"].astype(str).str.contains('CTNNB1', na = False, case = False), "APC"] = np.nan # in the APC column if any cell contains CTNNB1, replace with NaN
df_main.loc[df_main["CTNNB1"].astype(str).str.contains('APC', na = False, case = False), "CTNNB1"] = np.nan 

#######################
# COPY NUMBER DF
#######################
df_cn = df_main.copy()
df_cn = df_main.fillna(0)

# dictionary to map 
di = {"CN loss (deep deletion)": -2, "CN loss": -1, "CN loss (monoallelic)": -1, "Gain": 1, "Amp (low)": 2, "Amp (mid)": 2, "Amp (high)": 2}
df_cn = df_cn.replace(di)

# sometimes a sample has both CN alterations and mutations. this one accounts for that: 
for col in df_cn: df_cn.loc[df_cn[col].astype(str).str.contains('CN loss', na = False, case = False), col] = -1 # replace any string that contains "CN loss" with "-1, indicating monoallellic deletion
for col in df_cn: df_cn.loc[df_cn[col].astype(str).str.contains('Amplification', na = False, case = False), col] = 2

# drop the col with AR mutations 
df_cn = df_cn.drop("AR_muts", axis = 1)

# replace all strings with 0 because they have a CN of 0 
cols = df_cn.iloc[:, 4:].columns # get the cols we are interested in
mask = df_cn[cols].applymap(lambda x: isinstance(x, (int, float))) # true when float or integer
df_cn[cols] = df_cn[cols].where(mask).fillna(0) # we get na when we have a string, turn it into 0 

#######################
# MUTATION INFO DF
#######################
df_muts = df_main.copy()

# drop AR cn col and do some renaming 
df_muts = df_muts.drop("AR", axis = 1)
df_muts = df_muts.rename(columns = {"AR_muts" : "AR"})

# replace CN information with NAs
df_muts = df_muts.mask((df_muts == "CN loss (deep deletion)") | (df_muts == "PIK3CA amplification") | (df_muts == "APC CN loss (monoallelic)") | (df_muts == "CN loss") | (df_muts == "CN loss (monoallelic)") | (df_muts == "Gain") | (df_muts == "Amp (low)") | (df_muts == "Amp (mid)") | (df_muts == "Amp (high)"))

# sometimes a sample has both CN alterations and mutations. this one accounts for that: 
dict_to_replace = {"CN loss (monoallelic) and ": " ", " and CN loss": " ", "CN loss ": " ", "PIK3R1 CN loss": np.nan, 
                   "CN loss (deep deletion)": " ", "Gain": " ", "Amp (low)": " ", "Amp (mid)": " ", "Amp (high)": " "}

df_muts = df_muts.replace({" and CN loss": " ", "CN loss ": " ", "PIK3R1 CN loss": np.nan, "APC CN loss": np.nan}, regex = True)

searchfor = ["cn loss", "and amplification"]
for col in df_muts: df_muts.loc[df_muts[col].astype(str).str.contains('|'.join(searchfor), na = False, case = False), col] = np.nan # replace any string that contains "CN loss" with "-1, indicating monoallellic deletion

df_muts = df_muts.drop(["ctDNA fraction", "Mutation count", "DNA repair defect"], axis = 1) # drop unnecessary cols 
# for col in df_muts: df_muts.loc[df_muts[col].astype(str).str.contains(r'^(?=.*germline)(?=.*missense)', na = False, case = False), col] = np.nan # drop germline missense mutations 

# wide to long format
df_muts = df_muts.reset_index(drop = True)
df_muts = pd.melt(df_muts, id_vars='Sample', value_vars=df_muts.columns[df_muts.columns != 'Sample'].tolist())
df_muts = df_muts.dropna()

# sort based on sample id, reindex and rename the cols
df_muts = df_muts.sort_values(by=['Sample'])
df_muts = df_muts.reset_index(drop = True)
nameslist = ["Sample", "Gene", "Effect"]
df_muts.columns = nameslist

df_muts.loc[94, "Effect"] = "Intragenic arrangement"

# DEALING WITH ROWS WITH MULTIPLE EFFECTS IN ONE CELL
#######################

# some cells have more than one mutation in them, locate them by searching for "and": 
idx = df_muts.index[df_muts["Effect"].astype(str).str.contains("and", na = False, case = False)].tolist()
subsetted = df_muts.loc[idx, "Effect"].str.split("and ", expand = True) # choose muts repeated at least twice, split on "and" and expand into separate cols
subsetted = pd.concat([df_muts.loc[idx, ["Sample", "Gene"]], subsetted], axis=1) # concat to add sample name and gene name 
subsetted = pd.melt(subsetted, id_vars = ["Sample", "Gene"], value_vars = subsetted.iloc[:, 2:]).dropna()
del subsetted["vardf_cn.to_csv(, index = False)
iable"]
subsetted.columns = ["Sample", "Gene", "Effect"]

df_muts = df_muts.drop(idx, axis = 0) # drop the rows with "and" in them before we append the repeated version 
df_muts = pd.concat([df_muts, subsetted]).reset_index(drop = True)

# GERMLINE MUTS 
#######################
idx = df_muts.index[df_muts["Effect"].astype(str).str.contains("germline", na = False, case = False)].tolist() # find idx of germline muts 
df_germ = df_muts.loc[idx, :]
df_muts = df_muts.drop(idx, axis = 0) # remove the germline muts, will add them later
df_muts = df_muts.reset_index(drop = True)

df_germ["Effect"] = df_germ["Effect"].str.split(r" \(p").str[0] # strsplit to remove the mut location, and only keep the effect

# MUTS WITH GENE NAMES
#######################
# clean up mutations that have gene names in front of them 
df_muts = cleanup_gene("APC", df_muts, "Effect")
df_muts = cleanup_gene("CTNNB1", df_muts, "Effect")
df_muts = cleanup_gene("PIK3CA", df_muts, "Effect")

# REPEATED MUTS
#######################
# repeated muts have x2 in the strings
idx = df_muts.index[df_muts["Effect"].astype(str).str.contains("x2", na = False, case = False)].tolist() # find idx of repeated muts 
repeated = pd.concat([df_muts.iloc[idx, :]]*2) # repeat twice
repeated["Effect"] = repeated["Effect"].str.split(" x2").str[0] # remove the x2 string 
df_muts = df_muts.drop(idx, axis = 0) # drop the old row
df_muts = pd.concat([df_muts, repeated]) # add the repeated df

# FRAMESHIFT MUTS 
#######################colname
# search for strings that start with "(p" and replace that with "frameshift"
df_muts.loc[df_muts["Effect"].astype(str).str.startswith("p."), "Effect"] = "Frameshift mutation"
df_muts = remove_locations("Non-frameshift indel", "Non-frameshift indel", df_muts, "Effect")
df_muts = remove_locations("frameshift mutation", "Frameshift mutation", df_muts, "Effect")

# MISSENSE MUTS 
#######################
df_muts = remove_locations("missense", "Missense mutation", df_muts, "Effect")

# STOPGAIN MUTS 
#######################
df_muts = remove_locations("stopgain", "Stopgain mutation", df_muts, "Effect")
df_muts = remove_locations("second stopgain", "Second stopgain mutation", df_muts, "Effect")

# MULTIPLE SOMATIC MUTATIONS
df_muts["Effect"] == "Multiple somatic mutations"
idx = np.where(df_muts.loc[df_muts.Sample == "CAMP-0003_Baseline"]["Effect"] == "Multiple somatic mutations")

idx = df_muts["Effect"] == "Multiple somatic mutations"

df_muts = df_muts.drop(idx, )

df_muts.loc[df_muts.index.repeat(test.times)]



#######################
df_muts = pd.concat([df_muts, df_germ]) # add the germ line mutations we removed earlier 
df_muts['Effect'] = df_muts['Effect'].str.capitalize() # capitalize first letter of each mut effect, just in case
df_muts = df_muts.reset_index(drop = True) # reset index after adding the germline muts

df_muts = df_muts.replace("Mutliple somatic mutations", "Multiple somatic mutations")

# save the dfs as csv 
df_cn.to_csv("/groups/wyattgrp/users/amunzur/ind232/data/ind232_CN.csv", index = False)
df_muts.to_csv("/groups/wyattgrp/users/amunzur/ind232/data/ind232_muts.csv", index = False)