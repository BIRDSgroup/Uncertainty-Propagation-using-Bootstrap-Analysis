# -*- coding: utf-8 -*-
"""
Created on Sat Nov  8 20:47:57 2025

@author: HP
"""

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

parameter_names = ['$\\mathtt{obs}$', '$\\mu$', '$\\mu-\\sigma$', '$\\mu-2\\sigma$']

# %% The function
def fun(parameters: 'dict', reference: 'dict', klim: 'int', s: 'int', n: 'int', plottype: 'str' = 'cat', analysistype: 'str' = 'V', path: 'str' = None):
    x = pd.DataFrame(columns=parameter_names)
    for i in range(len(parameters)):
        j = 0 if analysistype == 'V' else i
        values = pd.Series(index=range(5,klim))        

        if plottype == 'recall':
            reference_set = reference[j].nlargest(klim).index

        for k in range(5,klim):
            if plottype == 'cat':
                reference_set = reference[j].nlargest(k).index
            parameter_set = parameters[i].nlargest(k).index

            values[k] = len(set(reference_set) & set(parameter_set)) / len(reference_set)
        x[parameter_names[i]] = values.copy()
    return x

# %% Main
centrality = 'pagerank'
klim = 1000
Tissues = ["Muscle_Skeletal", "Whole_Blood", "Skin_Sun_Exposed_Lower_leg", "Thyroid", "Lung"]

for s in [237, 73]:
    POG_vals = dict()
    for tissue in Tissues:
        print(tissue)
        order = pd.read_csv('../Gene Ordering/' + tissue + '_order.csv', header=None, index_col=None)[0]
        deg_pop_orig = pd.read_csv('../Centrality Computation/' + tissue + '/original_' + centrality + '.csv', header=None, index_col=None).iloc[0,:]
        n = len(deg_pop_orig)

        deg_obs = pd.read_csv('../Centrality Computation/' + tissue + '_' + str(s) + '/original_' + centrality + '.csv', header=None, index_col=None).iloc[0,:]

        deg_sample = pd.read_csv('../Centrality Computation/' + tissue + '_' + str(s) + '/sample_' + centrality + '.csv', index_col=0, header=None)
        deg_sample.columns = range(deg_sample.shape[1])

        deg_bootstrap = np.mean(deg_sample, axis = 0)
        std_bootstrap = np.std(deg_sample, axis = 0)
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
        POG_vals[tissue] = fun(parameters, parameters_pop, klim, s, n = len(order), plottype='cat', analysistype='V', path='')
    
    filename = 'Semi-Simulated Data - 1000/POG Values/' + centrality + '_' + str(s) + '.xlsx'
    with pd.ExcelWriter(filename, engine="openpyxl") as writer:
        for sheet, data in POG_vals.items():
            data.to_excel(writer, sheet_name=sheet, index=True, header=True)