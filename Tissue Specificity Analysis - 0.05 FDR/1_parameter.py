# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 18:12:10 2024

@author: HP

PARAMETER ANALYSIS
"""

import pandas as pd
import matplotlib.pyplot as plt
import os

# %% The function
def fun(centrality, Tissues, gtype):
    [os.makedirs('Parameter Analysis/' + gtype + '_' + centrality + '/' + tissue) for tissue in Tissues]
    names = ['deg', 'mu', 'mu-2sigma', 'mu-sigma']
    
    for tissue in Tissues:
        print(tissue)
        folder = tissue
        if '_' in tissue and tissue.split(sep='_')[-1].isdigit():
            folder = "_".join(tissue.split(sep='_')[:-1])
        
        all_genes = pd.read_csv('../Original Dataset/Preprocessed Files/' + folder + '/genes.csv', header=0, index_col=None, sep=',')
        specific_genes = pd.read_csv(gtype + ' Genes/' + folder + '.tsv', sep = '\t')
        indices = all_genes[all_genes['gene_name'].isin(specific_genes['Gene'])].index
        
        deg = pd.read_csv('../Centrality Computation/' + tissue + '/original_' + centrality + '.csv', header=None).iloc[0,:]        
        df = pd.read_csv('../Centrality Computation/' + tissue + '/sample_' + centrality + '.csv', header=None, index_col=0)
        df.columns = all_genes.index
        
        mu = df.mean(axis = 0)
        std = df.std(axis = 0)
        mu_2sigma = mu - (2 * std)
        mu_sigma = mu - std
        
        ranks = {
            0: deg.rank(method='first', ascending=False).astype('int'),
            1: mu.rank(method='first', ascending=False).astype('int'),
            2: mu_2sigma.rank(method='first', ascending=False).astype('int'),
            3: mu_sigma.rank(method='first', ascending=False).astype('int')
        }
        
        
        for i in range(1, 4):
            plt.figure()
            plt.scatter(ranks[0].drop(indices), ranks[i].drop(indices))
            plt.scatter(ranks[0][indices], ranks[i][indices])
            plt.axline((0, 0), slope=1, color='red')
            plt.xlabel('Ranking for ' + names[0])
            plt.ylabel('Ranking for ' + names[i])
            plt.title(names[0] + ' v/s '+ names[i] + '\n' + tissue)
            plt.savefig('Parameter Analysis/' + gtype + '_' + centrality + '/' + tissue + '/bias_' + names[i] + '.pdf')
            plt.close()
        
        for i in range(2, 4):
            plt.figure()
            plt.scatter(ranks[1].drop(indices), ranks[i].drop(indices))
            plt.scatter(ranks[1][indices], ranks[i][indices])
            plt.axline((0, 0), slope=1, color='red')
            plt.xlabel('Ranking for ' + names[1])
            plt.ylabel('Ranking for ' + names[i])
            plt.title(names[1] + ' v/s '+ names[i] + '\n' + tissue)
            plt.savefig('Parameter Analysis/' + gtype + '_' + centrality + '/' + tissue + '/mu vs ' + names[i] + '.pdf')
            plt.close()

# %% Main Function
Tissues = ['Kidney_Cortex', 'Lung', 'Muscle_Skeletal', 'Pancreas', 'Pituitary', 'Stomach', 'Thyroid', 'Whole_Blood', 'Muscle_Skeletal_73', 'Muscle_Skeletal_237', 'Whole_Blood_73', 'Whole_Blood_237', 'Lung_237', 'Lung_73', 'Vagina']

for gtype in ['Elevated', 'Enriched']:
    for centrality in ['degree', 'pagerank']:
        print(gtype, centrality)
        if gtype == 'Enriched':
            fun(centrality, Tissues[:-1], gtype)
        else:
            fun(centrality, Tissues, gtype)
