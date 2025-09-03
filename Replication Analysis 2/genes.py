# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 02:26:00 2025

@author: HP
"""

import numpy as np
import pandas as pd
import random

import qtl.io
import qtl.norm

import os
from time import time

tissue = 'Muscle_Skeletal'

# %% Recount Dataset
in_file = 'Original Dataset/' + tissue + '_recount/' + tissue + ' original.csv'
df_recount = pd.read_csv(in_file, sep = ',', header=0, index_col=0)
df_recount = df_recount.loc[df_recount['HG19en82 Gene Biotype'] == 'protein_coding', :]

# %% Normalization of the recount dataset
# TMM Normalization + Inverse Normal Transformation
df_recount_norm = df_recount.iloc[:,:-10]
n, s = df_recount_norm.shape
mask = (
        (np.sum(df_recount_norm >= 6, axis=1) >= 0.2 * s)
    ).values

tmm_counts_df = qtl.norm.edger_cpm(df_recount_norm, normalized_lib_sizes=True)
df_recount_norm = qtl.norm.inverse_normal_transform(tmm_counts_df[mask])

# %% GTEx Dataset
in_file = '../Original Dataset/Preprocessed Files/' + tissue + '/' + tissue + '.csv'
df_gtex = pd.read_csv(in_file, sep=',', header=0, index_col=0)

# %% Original Gene Lists
genes_recount = df_recount.loc[df_recount_norm.index, :]
genes_recount = genes_recount.iloc[:,-10:]
genes_recount = genes_recount.reset_index()
genes_gtex = pd.read_csv('../Original Dataset/Preprocessed Files/' + tissue + '/genes.csv')
genes_gtex['gene_id_new'] = genes_gtex['gene_id'].str.split('.').str[0]

df_gtex.set_index(genes_gtex['gene_id_new'], inplace=True)

# %% New Gene Lists
new_recount = genes_recount[ genes_recount['Tracking_ID'].isin(genes_gtex['gene_id_new']) ]
new_gtex = genes_gtex[ genes_gtex['gene_id_new'].isin(genes_recount['Tracking_ID']) ]

new_recount.set_index('Tracking_ID', inplace=True)
new_gtex.set_index('gene_id_new', inplace=True)

# %% Seed for ordering
seed_file = 'order_seed.txt'
if not os.path.exists(seed_file):
    seed = int(time())
    with open(seed_file, 'w') as file:
        print(seed, file = file)
    file.close()
else:
    with open(seed_file, 'r') as file:
        seed = int(file.read())
random.seed(seed)
np.random.seed(seed)

order = random.sample(new_gtex.index.tolist(), len(new_gtex))

# %% Extract in the same order
# GTEx
new_gtex_df = df_gtex.loc[order,:]
new_gtex = new_gtex.loc[order,:]

# Recount
new_recount_df = df_recount_norm.copy()
new_recount_df = new_recount_df.loc[order,:]
new_recount = new_recount.loc[order,:]

# %% Write to File
with open('gene ordering.csv','w') as order_file:
    print(*order, sep='\n', end='\n', file=order_file)
order_file.close()

# Recount Dataset
folder = 'Original Dataset/' + tissue + '_recount'
new_recount_df.to_csv(folder + '/' + tissue + '_recount normalized.csv', header=True, index=True)
new_recount.to_csv(folder + '/genes.csv', header=True, index=True)

# GTEx Dataset
folder = 'Original Dataset/' + tissue + '_gtex'
new_gtex_df.to_csv(folder + '/' + tissue + '_gtex.csv', header=True, index=True)
new_gtex.to_csv(folder + '/genes.csv', header=True, index=True)