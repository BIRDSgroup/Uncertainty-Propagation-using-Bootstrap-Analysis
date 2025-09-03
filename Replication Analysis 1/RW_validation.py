# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 14:45:55 2024

@author: HP
"""

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import sys

import plots

# %% Real - World Data
arg = sys.argv
tissue = arg[1]
centrality = arg[2]
s = int(arg[3])
klim = 300

print(s, tissue, sep = '\t')

# Read gene_ordering
order = pd.read_csv('../Gene Ordering/' + tissue + '_order.csv', header=None, index_col=None)[0]

# %% Read Centrality Metrics
deg_pop_orig = pd.read_csv('../Centrality Computation/' + tissue + '/original_' + centrality + '.csv', header=None, index_col=None).iloc[0,:]
n = len(deg_pop_orig)

deg_obs = pd.read_csv('../Centrality Computation/' + tissue + '_' + str(s) + '/original_' + centrality + '.csv', header=None, index_col=None).iloc[0,:]

deg_sample = pd.read_csv('../Centrality Computation/' + tissue + '_' + str(s) + '/sample_' + centrality + '.csv', index_col=0, header=None)
deg_sample.columns = range(deg_sample.shape[1])

deg_bootstrap = np.mean(deg_sample, axis = 0)
std_bootstrap = np.std(deg_sample, axis = 0)


# %%% Plots
print('Parameters')
parameters_pop = {
    0: deg_pop_orig[order],
    1: None,
    2: None,
    3: None
}

parameters = {
    0: deg_obs[order],
    1: deg_bootstrap[order],
    2: (deg_bootstrap - std_bootstrap)[order],
    3: (deg_bootstrap - (2 * std_bootstrap))[order]
}

# mu v/s sigma plot
filepath = 'Semi-Simulated Data/' + tissue + '/' + centrality + '_mu-sigma_' + str(s) + '.pdf'
plots.muVsSigma(deg_bootstrap[order], std_bootstrap[order], s, filepath)

# CAT Plot
filepath = 'Semi-Simulated Data/' + tissue + '/' + centrality + '_cat_' + str(s) + '.pdf'
plots.Plot(parameters, parameters_pop, klim, s, n = len(order), plottype='cat', analysistype='V', path=filepath)

# Recall @ k Plot
filepath = 'Semi-Simulated Data/' + tissue + '/' + centrality + '_recall_' + str(s) + '.pdf'
plots.Plot(parameters, parameters_pop, klim, s, n = len(order), plottype='recall', analysistype='V', path=filepath)
