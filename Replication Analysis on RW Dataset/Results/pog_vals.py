# -*- coding: utf-8 -*-
"""
Created on Sat Nov  8 20:47:57 2025

@author: HP
"""

import warnings
warnings.filterwarnings("ignore")

import os
import pandas as pd

parameter_names = ['$\\mathtt{obs}$', '$\\mu$', '$\\mu-\\sigma$', '$\\mu-2\\sigma$']

# %% The function
def fun(parameters: 'dict', reference: 'dict', klim: 'int', s: 'int', n: 'int', plottype: 'str' = 'cat', analysistype: 'str' = 'V', path: 'str' = None):
    x = pd.DataFrame(columns=parameter_names)
    for i in range(len(parameters)):
        j = 0 if analysistype == 'V' else i
        values = pd.Series(index=range(5,klim))        

        if plottype == 'recall':
            reference_set = reference[j].nlargest(klim).index

        for k in range(5,klim):
            if plottype == 'cat':
                reference_set = reference[j].nlargest(k).index
            parameter_set = parameters[i].nlargest(k).index

            values[k] = len(set(reference_set) & set(parameter_set)) / len(reference_set)
        x[parameter_names[i]] = values.copy()
    return x

# %% Main
centrality = 'pagerank'
klim = 1000
tissue = 'Muscle_Skeletal'

genes_gtex = pd.read_csv('../Original Dataset/Preprocessed Files/' + tissue + '_gtex/genes.csv')
genes_recount = pd.read_csv('../Original Dataset/Preprocessed Files/' + tissue + '_recount/genes.csv')
order = pd.read_csv('gene_order.csv', header=None)[0]

# GTEx
print('GTEx Dataset')
deg_gtex = pd.read_csv('../Centrality Computation/' + tissue + '_gtex/original_' + centrality + '.csv', header=None).iloc[0,:]
df_gtex = pd.read_csv('../Centrality Computation/' + tissue + '_gtex/sample_' + centrality + '.csv', header=None, index_col=0)
df_gtex.columns = deg_gtex.index

deg_gtex = deg_gtex[order]
df_gtex = df_gtex.loc[:, order]

mu_gtex = df_gtex.mean(axis = 0)
std_gtex = df_gtex.std(axis = 0)

parameters_gtex = {
    0: deg_gtex,
    1: mu_gtex,
    2: mu_gtex - std_gtex,
    3: mu_gtex - (2 * std_gtex)
}

# Recount
print('Recount Dataset')
deg_recount = pd.read_csv('../Centrality Computation/' + tissue + '_recount/original_' + centrality + '.csv', header=None).iloc[0,:]
df_recount = pd.read_csv('../Centrality Computation/' + tissue + '_recount/sample_' + centrality + '.csv', header=None, index_col=0)
df_recount.columns = deg_recount.index

deg_recount = deg_recount[order]
df_recount = df_recount.loc[:, order]

mu_recount = df_recount.mean(axis = 0)
std_recount = df_recount.std(axis = 0)

parameters_recount = {
    0: deg_recount,
    1: mu_recount,
    2: mu_recount - std_recount,
    3: mu_recount - (2 * std_recount)
}

POG_vals = fun(parameters_recount, parameters_gtex, klim, s = 53, n = len(order), plottype='cat', analysistype='R', path='')

filename = '1000/POG Values.xlsx'
mode = 'a' if os.path.exists(filename) else 'w'
with pd.ExcelWriter(filename, engine="openpyxl", mode=mode) as writer:
    POG_vals.to_excel(writer, sheet_name=centrality, index=True, header=True)