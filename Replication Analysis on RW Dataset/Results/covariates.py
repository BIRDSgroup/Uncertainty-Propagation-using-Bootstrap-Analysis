# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 18:51:08 2025

@author: HP
"""

import pandas as pd
import matplotlib.pyplot as plt

tissue = 'Muscle_Skeletal'
klim = 300

parameter_names = ['$\\mathtt{obs}$', '$\\mu$', '$\\mu-\\sigma$', '$\\mu-2\\sigma$']
file_names = ['obs', 'mu', 'mu-sigma', 'mu-2sigma']

# %% Main Part
for centrality in ['degree', 'pagerank']:
    # %% Random Order
    order = pd.read_csv('../../Replication Analysis 2 -- Wrong Results/Results/gene_order.csv', header=None)[0]
    
    # %% Without Covariate Adjustment -- Set 1
    deg_recount = pd.read_csv('../../Replication Analysis 2 -- Wrong Results/Centrality Computation/' + tissue + '_recount/original_' + centrality + '.csv', header=None).iloc[0,:]
    df_recount = pd.read_csv('../../Replication Analysis 2 -- Wrong Results/Centrality Computation/' + tissue + '_recount/sample_' + centrality + '.csv', header=None, index_col=0)
    df_recount.columns = deg_recount.index

    deg_recount = deg_recount[order]
    df_recount = df_recount.loc[:, order]

    mu_recount = df_recount.mean(axis = 0)
    std_recount = df_recount.std(axis = 0)

    parameters_recount_1 = {
        0: deg_recount,
        1: mu_recount,
        2: mu_recount - std_recount,
        3: mu_recount - (2 * std_recount)
    }
    
    genes_gtex = pd.read_csv('../../Replication Analysis 2 -- Wrong Results/Original Dataset/Preprocessed Files/' + tissue + '_gtex/genes.csv')
    genes_recount = pd.read_csv('../../Replication Analysis 2 -- Wrong Results/Original Dataset/Preprocessed Files/' + tissue + '_recount/genes.csv')    
    order = genes_gtex.loc[order, 'gene_id_new']
    
    # %% With Covariate Adjustment -- Set 2
    genes_gtex = pd.read_csv('../Original Dataset/Preprocessed Files/' + tissue + '_gtex/genes.csv')
    genes_recount = pd.read_csv('../Original Dataset/Preprocessed Files/' + tissue + '_recount/genes.csv')    
    order = pd.Index(genes_gtex.loc[:, 'gene_id_new']).get_indexer(order)
    
    deg_recount = pd.read_csv('../Centrality Computation/' + tissue + '_recount/original_' + centrality + '.csv', header=None).iloc[0,:]
    df_recount = pd.read_csv('../Centrality Computation/' + tissue + '_recount/sample_' + centrality + '.csv', header=None, index_col=0)
    df_recount.columns = deg_recount.index

    deg_recount = deg_recount[order]
    df_recount = df_recount.loc[:, order]

    mu_recount = df_recount.mean(axis = 0)
    std_recount = df_recount.std(axis = 0)

    parameters_recount_2 = {
        0: deg_recount,
        1: mu_recount,
        2: mu_recount - std_recount,
        3: mu_recount - (2 * std_recount)
    }
    
    # %% The Plots
    for i in range(len(parameter_names)):
        plt.figure()
        plt.scatter(parameters_recount_1[i], parameters_recount_2[i], color='#CC79A7')
        plt.axline((0,0), slope=1, color='#000000')
        plt.xlabel('Without Covariate Adjustment')
        plt.ylabel('With Covariate Adjustment')
        plt.title(parameter_names[i])
        plt.savefig('covs role/' + file_names[i] + '.pdf')
        plt.close()