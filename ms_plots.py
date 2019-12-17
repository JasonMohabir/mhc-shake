import numpy as np
import pandas as pd
from pathlib import Path
from mhcflurry import Class1AffinityPredictor
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt

# Load data and models - change path as necessary
DATA = Path("./data")
ms_file = DATA / "abelin_peptides.mhcflurry_no_mass_spec.csv"
MODELS = Path("/Users/haemin/Library/Application Support/mhcflurry/4/1.4.0")

no_ms_model = MODELS / "models_class1_selected_no_mass_spec/models"
no_ms_predictor = Class1AffinityPredictor.load(no_ms_model)

# Read data file and initial cleanup
ms = pd.read_csv(ms_file)
ms = ms.rename(columns={"mhcflurry": "mhcflurry2"})
ms["mhcflurry4"] = np.full_like(ms.mhcflurry2, -1)
ms = ms.loc[ms.allele.isin(no_ms_predictor.supported_alleles), :]
ms["peptide_len"] = ms.peptide.str.len()
ms.loc[ms["peptide_len"] > 12, "peptide_len"] = 13


'''
    Generate Figure A
'''
# Compute PPV values
models = ["netmhc", "netmhcpan", "mhcflurry2", "mhcflurry4"]
alleles = ms.allele[ms.allele.isin(no_ms_predictor.supported_alleles)].unique()

# Compute predictions for mhcflurry4
for allele in tqdm(alleles):
    curr_allele = ms[ms["allele"] == allele]
    ms.loc[curr_allele.index, "mhcflurry4"] = no_ms_predictor.predict(allele=allele, 
                                                                      peptides=curr_allele.peptide.values)


# Declare df's to hold computed PPV values

# PPV values for each allele, computed for each model tested
ppv_by_allele = pd.DataFrame(np.full((len(alleles), len(models)), -1.0), columns=models)
ppv_by_allele.index = pd.Series(alleles, dtype=str)
ppv_by_allele.index = ppv_by_allele.index.set_names(['allele'])

# PPV values for each allele, segmented by peptide length
peptide_lengths = ms["peptide_len"].unique()
index = pd.MultiIndex.from_product(np.array([alleles, peptide_lengths]), 
                                   names=['alleles', 'peptide_len'])
ppv_by_nmer = pd.DataFrame(np.full((index.size, len(models)), -1.0), index=index, columns=models)

# Populate df's
for allele in tqdm(alleles):
    curr_allele = ms[ms["allele"] == allele].copy()
    num_hits_allele = len(curr_allele[curr_allele.hit == 1].index.values)
    
    for model in models:            
        curr_model = curr_allele[["hit", model]]
        top_preds_model = curr_model.sort_values(by=model).iloc[:num_hits_allele]
        ppv_by_allele.loc[allele, model] = np.count_nonzero(top_preds_model.hit) / num_hits_allele
        
        grouped = curr_allele[[model, "hit", "peptide_len"]].groupby("peptide_len")
        for peptide_len, group in grouped:
            num_hits_nmer = len(group[group.hit == 1].index.values)
            top_preds_nmer = group.sort_values(by=model).iloc[:num_hits_nmer]
            ppv_by_nmer.loc[(allele, peptide_len), model] = np.count_nonzero(top_preds_nmer.hit) / num_hits_nmer

# Generate figure A
name_map = {
    "netmhc": "NetMHC 4.0", 
    "netmhcpan": "NetMHCpan 3.0", 
    "mhcflurry2": "MHCflurry 1.2.0 (no MS)",
    "mhcflurry4": "MHCflurry 1.4.0 (no MS)",
}

allele_df = ppv_by_allele.reset_index()
allele_df = allele_df.rename(columns=name_map)

by_model = [allele_df[[model, 'allele']].copy() for model in name_map.values()]
for name, df in zip(name_map.values(), by_model):
    df['model'] = pd.Series([name]*len(df))
    df.rename({name: "ppv"}, axis=1, inplace=True)

ppv_df = pd.concat(by_model, axis=0, ignore_index=True)
ppv_df['allele_name'] = ppv_df['allele'].str.split('-').str[1]

ppv_df['pos'] = np.zeros((len(ppv_df.index.values), 1))
for i, name in enumerate(name_map.values()):
    ppv_df.loc[ppv_df.model == name, 'pos'] = i
ppv_df = ppv_df.sort_values(by=['pos', 'allele'])

fig, ax = plt.subplots(1,1, figsize=(15,10))
ax.grid(b=True, which='major')
plt.xticks(rotation=90)
sns.barplot(data=ppv_df, x="allele_name", y='ppv', hue="model", ax=ax, ci="sd", errcolor="k", errwidth=1)
plt.legend(loc='upper right')
plt.show()
plt.close()

# Generate figure B
ppv_by_nmer['mhcflurry2_netmhcpan'] = ((ppv_by_nmer['mhcflurry2'] / ppv_by_nmer['netmhcpan'])-1)*100
ppv_by_nmer['mhcflurry4_netmhcpan'] = ((ppv_by_nmer['mhcflurry4'] / ppv_by_nmer['netmhcpan'])-1)*100
ppv_by_nmer['mhcflurry4_mhcflurry2'] = ((ppv_by_nmer['mhcflurry4'] / ppv_by_nmer['mhcflurry2'])-1)*100
ppv_by_nmer.replace(to_replace=np.inf, value=np.nan, inplace=True)
ppv_by_nmer.dropna(inplace=True)

ppv_by_len = ppv_by_nmer.reset_index(level=['peptide_len', 'alleles'])
ppv_by_len.drop('alleles', axis=1, inplace=True)

all_lengths = ppv_by_len.copy().groupby('peptide_len').aggregate(np.median)
all_lengths = all_lengths.reset_index(level=['peptide_len'])
all_lengths['peptide_len'] = 6

non_9 = ppv_by_len[ppv_by_len["peptide_len"] != 9].copy().groupby('peptide_len').aggregate(np.mean)
non_9 = non_9.reset_index(level=['peptide_len'])
non_9['peptide_len'] = 7

nmer_df = pd.concat([ppv_by_len, all_lengths, non_9], axis=0, ignore_index=True)
nmer_df['nmer'] = nmer_df.peptide_len.apply(str)
nmer_df['nmer'] = nmer_df['nmer'].apply(lambda x: x+"-mers")
nmer_df.loc[nmer_df.peptide_len == 6, 'nmer'] = 'all lengths'
nmer_df.loc[nmer_df.peptide_len == 7, 'nmer']= 'non-9-mers'
nmer_df.loc[nmer_df.peptide_len == 13, 'nmer']= '>12-mers'
nmer_df.sort_values(by='peptide_len', inplace=True)

# Plot figure for v1.2 vs NetMHCpan
fig, ax = plt.subplots(1,1, figsize=(10,5))
ax.set_xlim(-40, 100)
ax.grid(b=True, which='major')
sns.boxplot(data=nmer_df, orient='h', x='mhcflurry2_netmhcpan', y='nmer', ax=ax, whis=0)
ax.set_ylabel("")
ax.set_xlabel("Improvement of MHCFlurry 1.2.0 over NetMHCpan 3.0\n(%PPV)")
plt.show()
plt.close()


# Plot figure for v1.4 vs NetMHCpan
fig, ax = plt.subplots(1,1, figsize=(10,5))
ax.set_xlim(-40, 100)
ax.grid(b=True, which='major')
sns.boxplot(data=nmer_df, orient='h', x='mhcflurry4_netmhcpan', y='nmer', ax=ax, whis=0)
ax.set_ylabel("")
ax.set_xlabel("Improvement of MHCFlurry 1.4.0 over NetMHCpan 3.0\n(%PPV)")
plt.show()
plt.close()

