# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 18:51:08 2025

@author: HP
"""

import pandas as pd

tissue = 'Muscle_Skeletal'
klim = 1000

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
        
    genes_gtex.set_index('gene_id_new', inplace=True)
    genes_recount.set_index('Tracking_ID', inplace=True)
    order = genes_gtex.index[order]

    # %% GTEx
    print('GTEx Dataset')
    deg_gtex = pd.read_csv('../Centrality Computation/' + tissue + '_gtex/original_' + centrality + '.csv', header=None).iloc[0,:]
    df_gtex = pd.read_csv('../Centrality Computation/' + tissue + '_gtex/sample_' + centrality + '.csv', header=None, index_col=0)
    
    deg_gtex.index = genes_gtex.index
    df_gtex.columns = genes_gtex.index

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
    for k in parameters_gtex.keys():
        parameters_gtex[k] = parameters_gtex[k].rank(method='first', ascending=False)

    # %% Recount
    print('Recount Dataset')
    deg_recount = pd.read_csv('../Centrality Computation/' + tissue + '_recount/original_' + centrality + '.csv', header=None).iloc[0,:]
    df_recount = pd.read_csv('../Centrality Computation/' + tissue + '_recount/sample_' + centrality + '.csv', header=None, index_col=0)
    
    deg_recount.index = genes_recount.index
    df_recount.columns = genes_recount.index
    
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
    for k in parameters_recount.keys():
        parameters_recount[k] = parameters_recount[k].rank(method='first', ascending=False)
    
    # %% The Tables
    print("The Tables")
    specific_genes = pd.read_csv(tissue + '.tsv', header=0, index_col=None, sep='\t')
    specific_genes.set_index('Ensembl', inplace=True)
    top_gene_list = specific_genes['RNA tissue specificity score'].nlargest(17).index
    top_gene_list = top_gene_list.intersection(genes_gtex.index)
    
    top_genes_table = pd.DataFrame({
        'Gene Name': specific_genes.loc[top_gene_list, 'Gene'],
        'Specificity Type': specific_genes.loc[top_gene_list, 'RNA tissue specificity'],
        'Specificity Score': specific_genes.loc[top_gene_list, 'RNA tissue specificity score']
    })
    for i in range(4):
        top_genes_table[str(i) + '_disc'] = parameters_gtex[i][top_gene_list]
        top_genes_table[str(i) + '_repl'] = parameters_recount[i][top_gene_list]
    
    top_genes_table.index.name = 'Ensembl ID'
    top_genes_table.to_csv(str(klim) + '/' + centrality + '.csv', index=True, header=True)
    