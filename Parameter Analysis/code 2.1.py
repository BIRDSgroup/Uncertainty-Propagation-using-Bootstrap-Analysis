# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 13:38:33 2025

@author: HP
"""

import sys
args = sys.argv
centrality = args[1]

# %% Packages
import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

Tissues = pd.read_excel('../Research Plan.xlsx', sheet_name='Sheet2', header=0)

names = ['obs', 'mu', 'mu-sigma', 'mu-2sigma']
colors = ['#000000', '#0072B2', '#E69F00', '#009E73']
shapes = ['o', '^', 's', '*']
plt.style.use("tableau-colorblind10")

# %% mu v/s mu-sigma v/s mu-2sigma

corr = pd.DataFrame(columns=['0-1', '1-2', '1-3', '2-3'], index=Tissues['Tissue Name'])
for i in Tissues.index:
    print(Tissues.loc[i,'Tissue Name'])
    deg = pd.read_csv('../Centrality Computation/' + Tissues.loc[i,'Folder Name'] + '/original_' + centrality + '.csv', header=None).iloc[0,:]
    df = pd.read_csv('../Centrality Computation/' + Tissues.loc[i,'Folder Name'] + '/sample_' + centrality + '.csv', header=None, index_col=0)
    df.columns = range(df.shape[1])
    
    mu = df.mean(axis = 0)
    std = df.std(axis = 0)
    mu_2sigma = mu - (2 * std)
    mu_sigma = mu - std
    
    parameters = {
        0: deg,
        1: mu,
        2: mu_sigma,
        3: mu_2sigma
    }
    
    corr.loc[Tissues.loc[i,'Tissue Name'], '0-1'] = spearmanr(parameters[0], parameters[1]).statistic
    corr.loc[Tissues.loc[i,'Tissue Name'], '1-2'] = spearmanr(parameters[1], parameters[2]).statistic
    corr.loc[Tissues.loc[i,'Tissue Name'], '1-3'] = spearmanr(parameters[1], parameters[3]).statistic
    corr.loc[Tissues.loc[i,'Tissue Name'], '2-3'] = spearmanr(parameters[2], parameters[3]).statistic
    
# %%% Plots
labels = ['Correlation between $\\mathtt{obs}$ and $\\mu$', 'Correlation between $\\mu$ and $\\mu-\\sigma$', 'Correlation between $\\mu$ and $\\mu-2\\sigma$', 'Correlation between $\\mu-\\sigma$ and $\\mu-2\\sigma$']
plt.figure()
for i in range(4):
    plt.scatter(Tissues['# Samples'], corr[corr.columns[i]], label = labels[i], marker=shapes[i], color=colors[i])
plt.legend(loc = 'upper left', bbox_to_anchor=(0.35, -0.25))
plt.xlabel('Sample Size ($s$)')
plt.ylabel('Correlation Coefficient ($r$)')
plt.title('Spearman Correlation among Centrality Metrics \n (with respect to sample size)')
plt.tight_layout()
plt.savefig('combined plots/' + centrality + '_all tissues.pdf')
plt.close()
