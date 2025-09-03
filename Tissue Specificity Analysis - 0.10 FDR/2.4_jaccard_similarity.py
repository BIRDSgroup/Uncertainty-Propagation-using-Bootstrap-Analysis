#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  3 15:49:45 2025

@author: sugyani
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# %% Main Function
Pathways = ['Muscle_Skeletal', 'Lung', 'Thyroid', 'Pancreas', 'Pituitary', 'Stomach', 'Kidney_Cortex', 'Vagina']

for gtype in ['Elevated', 'Enriched']:
    if gtype == 'Enriched':
        Pathways = Pathways[:-1]
    JC = pd.DataFrame(columns=Pathways, index=Pathways, dtype=float)
    for path1 in Pathways:
        genes1 = set(pd.read_csv(gtype + ' Genes/' + path1 + '.tsv', sep = '\t')['Gene'])
        for path2 in Pathways:
            genes2 = set(pd.read_csv(gtype + ' Genes/' + path2 + '.tsv', sep = '\t')['Gene'])
            JC.loc[path1, path2] = JC.loc[path2, path1] = len(genes1 & genes2) / len(genes1 | genes2)
    
    plt.figure()
    ax = sns.heatmap(JC, cmap='coolwarm', vmin=-1, vmax=1, xticklabels=True, yticklabels=True)
    ax.hlines([3, 6], *ax.get_xlim(), colors='white', linewidths=2)
    ax.vlines([3, 6], *ax.get_ylim(), colors='white', linewidths=2)
    plt.xticks(rotation=45, ha='right')
    plt.title('Similarity between Tissues')
    plt.tight_layout()
    plt.savefig('Pathway Results/' + gtype + '.pdf')
    plt.close()
    
    del JC