# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 18:12:10 2024

@author: HP

PARAMETER ANALYSIS
"""

import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

all_Tissues = pd.read_excel('../../Research Plan.xlsx', sheet_name='Sheet2', header=0, index_col=1)
parameter_names = ['$\\deg$', '$\\mu$', '$\\mu-\\sigma$', '$\\mu-2\\sigma$']
colors = ['#000000', '#0072B2', '#E69F00', '#009E73']
shapes = ['o', '^', 's', '*']
annot_col = ['#56B4E9', '#D55E00', '#CC79A7']
plt.style.use("tableau-colorblind10")

klim = 100

# %% The function
def fun(centrality, Tissues, gtype):
    for tissue in Tissues:
        print(tissue)
        folder = tissue
        if '_' in tissue and tissue.split(sep='_')[-1].isdigit():
            folder = "_".join(tissue.split(sep='_')[:-1])                    
        all_genes = pd.read_csv('../../Original Dataset/Preprocessed Files/' + folder + '/genes.csv', header=0, index_col=None, sep=',')
        specific_genes = pd.read_csv('../' + gtype + ' Genes/' + folder + '.tsv', sep = '\t')
        indices = all_genes[all_genes['gene_name'].isin(specific_genes['Gene'])].index
        
        order = pd.read_csv('../../Gene Ordering/' + folder + '_order.csv', header=None, index_col=None)[0]
        
        deg = pd.read_csv('../../Centrality Computation/' + tissue + '/original_' + centrality + '.csv', header=None).iloc[0,:]        
        df = pd.read_csv('../../Centrality Computation/' + tissue + '/sample_' + centrality + '.csv', header=None, index_col=0)
        df.columns = all_genes.index
        
        mu = df.mean(axis = 0)
        std = df.std(axis = 0)
        mu_2sigma = mu - (2 * std)
        mu_sigma = mu - std
        
        parameters = {
            0: deg[order],
            1: mu[order],
            2: mu_sigma[order],
            3: mu_2sigma[order]
        }
        
        # Recall @ k Plot
        plt.figure()
        for i in range(4):
            recall = pd.Series(index=range(5,klim))
            deg_pop_indices = np.argsort(deg_pop)[-klim:]
            for k in range(5,klim):
                parameter_i_indices = np.argsort(parameters[i])[-k:]
                recall[k] = len(set(deg_pop_indices) & set(parameter_i_indices))
            plt.plot(recall, label=parameter_names[i], color=colors[i])
        plt.legend()

        
    
# %% Main Function
Tissues = ['Kidney_Cortex', 'Lung', 'Muscle_Skeletal', 'Pancreas', 'Pituitary', 'Stomach', 'Thyroid', 'Whole_Blood', 'Muscle_Skeletal_73', 'Muscle_Skeletal_237', 'Whole_Blood_73', 'Whole_Blood_237', 'Lung_237', 'Lung_73', 'Vagina']

for gtype in ['Elevated', 'Enriched']:
    for centrality in ['degree', 'pagerank']:
        print(gtype, centrality)
        if gtype == 'Enriched':
            fun(centrality, Tissues[:-1], gtype)
        else:
            fun(centrality, Tissues, gtype)

