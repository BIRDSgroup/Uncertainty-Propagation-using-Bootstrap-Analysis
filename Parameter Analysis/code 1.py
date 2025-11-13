# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 12:57:42 2024

@author: HP
"""


import pandas as pd
import matplotlib.pyplot as plt
import os

Tissues = ['Brain_Amygdala', 'Brain_Caudate_basal_ganglia', 'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Substantia_nigra', 'Kidney_Cortex', 'Lung', 'Muscle_Skeletal', 'Pancreas', 'Pituitary', 'Skin_Sun_Exposed_Lower_leg', 'Stomach', 'Thyroid', 'Uterus', 'Vagina', 'Whole_Blood', 'Muscle_Skeletal_73', 'Muscle_Skeletal_237', 'Whole_Blood_73', 'Whole_Blood_237', 'Lung_237', 'Lung_73']
centrality = 'pagerank'

[os.makedirs('rank plots/' + centrality + '/' + tissue) for tissue in Tissues]
[os.makedirs('value plots/' + centrality + '/' + tissue) for tissue in Tissues]
[os.makedirs('ranks/' + centrality + '/' + tissue) for tissue in Tissues]

names = ['deg', 'mu', 'mu-2sigma', 'mu-sigma']

# %% mu v/s mu-sigma v/s mu-2sigma
for tissue in Tissues:
    print(tissue)
    deg = pd.read_csv('../Centrality Computation/' + tissue + '/original_' + centrality + '.csv', header=None).iloc[0,:]
    df = pd.read_csv('../Centrality Computation/' + tissue + '/sample_' + centrality + '.csv', header=None, index_col=0)
    df.columns = range(df.shape[1])
    
    mu = df.mean(axis = 0)
    std = df.std(axis = 0)
    mu_2sigma = mu - (2 * std)
    mu_sigma = mu - std
    
    parameters = {
        0: deg,
        1: mu,
        2: mu_2sigma,
        3: mu_sigma
    }
    
    ranks = {
        0: deg.rank(method='first', ascending=False).astype('int'),
        1: mu.rank(method='first', ascending=False).astype('int'),
        2: mu_2sigma.rank(method='first', ascending=False).astype('int'),
        3: mu_sigma.rank(method='first', ascending=False).astype('int')
    }
    
    
    for i in range(1, 4):
        plt.figure()
        plt.scatter(ranks[0], ranks[i])
        plt.axline((0, 0), slope=1, color='red')
        plt.xlabel('Ranking for ' + names[0])
        plt.ylabel('Ranking for ' + names[i])
        plt.title(names[0] + ' v/s '+ names[i] + '\n' + tissue)
        plt.savefig('rank plots/' + centrality + '/' + tissue + '/bias_' + names[i] + '.pdf')
        plt.close()

    for i in range(2, 4):
        plt.figure()
        plt.scatter(ranks[1], ranks[i])
        plt.axline((0, 0), slope=1, color='red')
        plt.xlabel('Ranking for ' + names[1])
        plt.ylabel('Ranking for ' + names[i])
        plt.title(names[1] + ' v/s '+ names[i] + '\n' + tissue)
        plt.savefig('rank plots/' + centrality + '/' + tissue + '/mu vs ' + names[i] + '.pdf')
        plt.close()
    
    for i in range(1, 4):
        plt.figure()
        plt.scatter(parameters[0], parameters[i])
        plt.axline((0, 0), slope=1, color='red')
        plt.xlabel('Value for ' + names[0])
        plt.ylabel('Value for ' + names[i])
        plt.title(names[0] + ' v/s '+ names[i] + '\n' + tissue)
        plt.savefig('value plots/' + centrality + '/' + tissue + '/bias_' + names[i] + '.pdf')
        plt.close()

    for i in range(2, 4):
        plt.figure()
        plt.scatter(parameters[1], parameters[i])
        plt.axline((0, 0), slope=1, color='red')
        plt.xlabel('Value for ' + names[1])
        plt.ylabel('Value for ' + names[i])
        plt.title(names[1] + ' v/s '+ names[i] + '\n' + tissue)
        plt.savefig('value plots/' + centrality + '/' + tissue + '/mu vs ' + names[i] + '.pdf')
        plt.close()
        
    folder = tissue
    if '_' in tissue and tissue.split(sep='_')[-1].isdigit():
        folder = "_".join(tissue.split(sep='_')[:-1])
    genes = pd.read_csv('../Original Dataset/Preprocessed Files/' + folder + '/genes.csv', header=0, index_col=None, sep=',')
    output = pd.concat([genes['gene_name'], deg, mu, mu_2sigma, mu_sigma], axis = 1)
    output.columns = ['genes'] + names
    for p in names:
        opfile = 'ranks/' + centrality + '/' + tissue + '/' + tissue + '_' + p + '.rnk'
        (output.loc[:, ['genes', p]]).to_csv(opfile, sep = '\t', header = False, index=False)
    
