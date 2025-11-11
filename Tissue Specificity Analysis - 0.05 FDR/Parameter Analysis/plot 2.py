# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 14:45:55 2024

@author: HP
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

names = ['$\\deg$', '$\\mu$', '$\\mu-\\sigma$', '$\\mu-2\\sigma$']
colors = ['#000000', '#0072B2', '#E69F00', '#009E73']
shapes = ['o', '^', 's', '*']
annot_col = ['#56B4E9', '#D55E00', '#CC79A7']

# %% The function
def fun(centrality, Tissues, gtype):
    r_pop = pd.DataFrame(columns=names + ['s'])
    r_obs = pd.DataFrame(columns=names + ['s'])
    for tissue in Tissues:
        print(tissue)
        all_genes = pd.read_csv('../../Original Dataset/Preprocessed Files/' + tissue + '/genes.csv', header=0, index_col=None, sep=',')
        specific_genes = pd.read_csv('../' + gtype + ' Genes/' + tissue + '.tsv', sep = '\t')
        indices = all_genes[all_genes['gene_name'].isin(specific_genes['Gene'])].index
        
        for s in [73, 237]:
            print(s)
            deg_pop = pd.read_csv('../../Centrality Computation/' + tissue + '/original_' + centrality + '.csv', header=None, index_col=None).iloc[0,:]            
            deg_obs = pd.read_csv('../../Centrality Computation/' + tissue + '_' + str(s) + '/original_' + centrality + '.csv', header=None, index_col=None).iloc[0,:]

            deg_sample = pd.read_csv('../../Centrality Computation/' + tissue + '_' + str(s) + '/sample_' + centrality + '.csv', index_col=0, header=None)
            deg_sample.columns = range(deg_sample.shape[1])

            deg_bootstrap = np.mean(deg_sample, axis = 0)
            std_bootstrap = np.std(deg_sample, axis = 0)
            
            parameters = {
                1: deg_bootstrap,
                2: deg_bootstrap - std_bootstrap,
                3: deg_bootstrap - (2 * std_bootstrap)
            }
            
            for i in range(1,4):
                r_pop.loc[tissue + '_' + str(s), [names[i], 's']] = [spearmanr(parameters[i][indices], deg_pop[indices]).statistic, s]
                r_obs.loc[tissue + '_' + str(s), [names[i], 's']] = [spearmanr(parameters[i][indices], deg_obs[indices]).statistic, s]
    
    # Population as Ground Truth Plot
    plt.figure()
    for i in range(1,4):
        plt.scatter(r_pop['s'], r_pop[names[i]], label = names[i], color = colors[i], marker=shapes[i], s = 10)
    plt.xlim(0, 300)
    plt.xticks([73, 237])
    plt.legend()
    plt.xlabel('Sample Size ($s$)')
    plt.ylabel('Correlation Coefficient ($r$)')
    for i, tissue in enumerate(Tissues):
        for s in [73, 237]:
            for parameter in names[1:4]:
                point = (s, r_pop.loc[tissue + '_' + str(s), parameter])
                plt.annotate(tissue, xy=point, xytext=(125, 0.35+i*0.1), color=annot_col[i], arrowprops=dict(arrowstyle="->", color=annot_col[i]))
    plt.title('Correlation of Tissue-Specific Genes \n (Population Degree as Ground Truth)')
    plt.savefig('Population_' + gtype + '_' + centrality + '.pdf')
    plt.close()
    
    # Observed as Ground Truth Plot
    plt.figure()
    for i in range(1,4):
        plt.scatter(r_obs['s'], r_obs[names[i]], label = names[i], color = colors[i], marker=shapes[i], s = 10)
    plt.xlim(0, 300)
    plt.xticks([73, 237])
    plt.legend()
    plt.xlabel('Sample Size ($s$)')
    plt.ylabel('Correlation Coefficient ($r$)')
    for i, tissue in enumerate(Tissues):
        for s in [73, 237]:
            for parameter in names[1:4]:
                point = (s, r_obs.loc[tissue + '_' + str(s), parameter])
                plt.annotate(tissue, xy=point, xytext=(125, 0.35+i*0.1), color=annot_col[i], arrowprops=dict(arrowstyle="->", color=annot_col[i]))
    plt.title('Correlation of Tissue-Specific Genes \n (Observed Degree as Ground Truth)')
    plt.savefig('Observed_' + gtype + '_' + centrality + '.pdf')
    plt.close()
    
# %% Main Function
Tissues = ['Muscle_Skeletal', 'Whole_Blood', 'Lung']

for gtype in ['Elevated', 'Enriched']:
    for centrality in ['degree', 'pagerank']:
        print(gtype, centrality)
        fun(centrality, Tissues, gtype)


