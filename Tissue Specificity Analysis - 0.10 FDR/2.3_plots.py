# -*- coding: utf-8 -*-
"""
Created on Sat Jan  4 19:52:00 2025

@author: HP
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

names = ['$\\deg$', '$\\mu$', '$\\mu-\\sigma$', '$\\mu-2\\sigma$']
colors = ['#7f7f7f', '#1f77b4', '#ff7f0e', '#2ca02c']

# %% Heatmaps
def plot_hmap(centrality, gtype, Tissues, Paths):
    parameters = ['deg', 'mu', 'mu_sigma', 'mu_2sigma']
    for i, parameter in enumerate(parameters):
        print(parameter)
        heatmap_data = pd.DataFrame(columns= Paths, index = Tissues, dtype=float) # rows are tissues and columns are paths
        
        for tissue in Tissues:
            filepath = 'Pathway Analysis/' + gtype + '_' + centrality + '/' + tissue + '/Project_' + parameter + '/enrichment_results_' + parameter + '.txt'
            if not os.path.exists(filepath):
                continue
            results = pd.read_csv(filepath, sep='\t')
            
            enriched_filepath = 'Pathway Analysis/' + gtype + '_' + centrality + '/' + tissue + '/Project_' + parameter + '/enriched_geneset_wsc_topsets_' + parameter + '.txt'
            enriched_paths = pd.read_csv(enriched_filepath).iloc[:,0]
            
            for path in enriched_paths:
                heatmap_data.loc[tissue, path] = results.loc[results['geneSet'] == path, 'enrichmentScore'].values
        
        heatmap_data.fillna(0, inplace=True)
        
        plt.figure()
        ax = sns.heatmap(heatmap_data, cmap='coolwarm', vmin=-1, vmax=1, xticklabels=True, yticklabels=True)
        ax.hlines([4, 9], *ax.get_xlim(), colors='white', linewidths=2)
        ax.vlines([4, 7], *ax.get_ylim(), colors='white', linewidths=2)
        plt.xticks(rotation=45, ha='right')
        plt.title('All Paths v/s All Tissues -- ' + names[i])
        plt.xlabel('Pathways')
        plt.ylabel('Tissues')
        plt.tight_layout()
        plt.savefig('Pathway Results/Heatmaps/' + gtype + '/' + centrality + '_' + parameter + '.pdf')
        plt.close()
    
# %% Main Function
Tissues = ['Whole_Blood', 'Muscle_Skeletal', 'Lung', 'Thyroid', 'Pancreas', 'Pituitary', 'Stomach', 'Muscle_Skeletal_237', 'Lung_237', 'Kidney_Cortex', 'Vagina', 'Muscle_Skeletal_73',  'Lung_73']
Paths = ['Muscle_Skeletal', 'Lung', 'Thyroid', 'Pancreas', 'Pituitary', 'Stomach', 'Kidney_Cortex', 'Vagina']

for gtype in ['Elevated', 'Enriched']:
    os.makedirs('Pathway Results/Heatmaps/' + gtype)
    os.makedirs('Pathway Results/Recall/' + gtype)
    for centrality in ['degree', 'pagerank']:
        print(gtype, centrality)
        plot_hmap(centrality, gtype, Tissues, Paths)