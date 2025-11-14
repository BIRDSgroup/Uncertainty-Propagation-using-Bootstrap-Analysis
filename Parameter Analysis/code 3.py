# -*- coding: utf-8 -*-
"""
Created on Sat Nov 15 02:57:24 2025

@author: HP
"""

import pandas as pd
import os

Tissues = ['Brain_Amygdala', 'Brain_Caudate_basal_ganglia', 'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Substantia_nigra', 'Kidney_Cortex', 'Lung', 'Muscle_Skeletal', 'Pancreas', 'Pituitary', 'Skin_Sun_Exposed_Lower_leg', 'Stomach', 'Thyroid', 'Uterus', 'Vagina', 'Whole_Blood', 'Muscle_Skeletal_73', 'Muscle_Skeletal_237', 'Whole_Blood_73', 'Whole_Blood_237', 'Skin_Sun_Exposed_Lower_leg_237', 'Skin_Sun_Exposed_Lower_leg_73', 'Thyroid_237', 'Thyroid_73', 'Lung_237', 'Lung_73']
centrality = 'pagerank'

[os.makedirs('metric values/' + centrality + '/' + tissue) for tissue in Tissues]

names = ['obs', 'mu', 'mu-sigma', 'mu-2sigma', 'std']

for tissue in Tissues:
    print(tissue)
    deg = pd.read_csv('../Centrality Computation/' + tissue + '/original_' + centrality + '.csv', header=None).iloc[0,:]
    df = pd.read_csv('../Centrality Computation/' + tissue + '/sample_' + centrality + '.csv', header=None, index_col=0)
    df.columns = range(df.shape[1])
    
    folder = tissue
    if '_' in tissue and tissue.split(sep='_')[-1].isdigit():
        folder = "_".join(tissue.split(sep='_')[:-1])
    genes = pd.read_csv('../Original Dataset/Preprocessed Files/' + folder + '/genes.csv', header=0, index_col=None, sep=',')
    
    mu = df.mean(axis = 0)
    std = df.std(axis = 0)
    
    parameters = {
        0: deg,
        1: mu,
        2: (mu - std),
        3: (mu - (2*std)),
        4: std
    }
    
    op_folder = 'metric values/' + centrality + '/' + tissue + '/'
    for i in range(len(parameters)):
        parameters[i].index = genes['gene_name']
        filepath = op_folder + names[i] + '.csv'
        parameters[i].to_csv(filepath, header=False, index=True, sep=',')
