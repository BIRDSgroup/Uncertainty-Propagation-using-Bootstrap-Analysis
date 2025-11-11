# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 18:51:08 2025

@author: HP
"""

import pandas as pd

import plots


tissue = 'Muscle_Skeletal'
klim = 300

# %% Main Part
for centrality in ['degree', 'pagerank']:
    # %% Random Order
    import os
    import random
    from time import time
    
    genes_gtex = pd.read_csv('../Original Dataset/Preprocessed Files/' + tissue + '_gtex/genes.csv')
    genes_recount = pd.read_csv('../Original Dataset/Preprocessed Files/' + tissue + '_recount/genes.csv')    
    if not os.path.exists('gene_order.csv'):
        seed_file = 'seed.txt'
        if not os.path.exists(seed_file):
            seed = int(time())
            with open(seed_file, 'w') as file:
                print(seed, file = file)
            file.close()
        else:
            with open(seed_file, 'r') as file:
                seed = int(file.read())
        random.seed(seed)
        #np.random.seed(seed)

        order = random.sample(range(len(genes_gtex)), len(genes_gtex))
        with open('gene_order.csv', 'w') as file:
            print(*order, sep = '\n', end='\n', file=file)
        file.close()
    else:
        order = pd.read_csv('gene_order.csv', header=None)[0]

    # %% GTEx
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

    # %% Recount
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


    # %% Plots
    print('Plotting')
    # CAT Plot
    filepath = str(klim) + '/' + centrality + '_cat.pdf'
    plots.Plot(parameters_recount, parameters_gtex, klim, s = 53, n = len(order), plottype='cat', analysistype='R', path=filepath)

    # Recall Plot
    filepath = str(klim) + '/' + centrality + '_recall.pdf'
    plots.Plot(parameters_recount, parameters_gtex, klim, s = 53, n = len(order), plottype='recall', analysistype='R', path=filepath)