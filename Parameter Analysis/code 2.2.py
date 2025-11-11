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

Tissues = pd.read_excel('../Research Plan.xlsx', sheet_name='Sheet2', header=0)

# %% Subsampled Tissues
labels = ['Correlation between $\\mathtt{obs}$ and $\\mu$', 'Correlation between $\\mu$ and $\\mu-\\sigma$', 'Correlation between $\\mu$ and $\\mu-2\\sigma$', 'Correlation between $\\mu-\\sigma$ and $\\mu-2\\sigma$']
corr = pd.DataFrame(columns=['s', '0-1', '1-2', '1-3', '2-3'], index=['Muscle_Skeletal_237', 'Whole_Blood_237', 'Skin_Sun_Exposed_Lower_leg_237', 'Thyroid_237', 'Lung_237', 'Muscle_Skeletal_73', 'Whole_Blood_73', 'Skin_Sun_Exposed_Lower_leg_73', 'Thyroid_73', 'Lung_73'])

for tissue in corr.index:
    print(tissue)
    deg = pd.read_csv('../Centrality Computation/' + tissue + '/original_' + centrality + '.csv', header=None).iloc[0,:]
    df = pd.read_csv('../Centrality Computation/' + tissue + '/sample_' + centrality + '.csv', header=None, index_col=0)
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
    
    corr.loc[tissue, 's'] = int(tissue.split('_')[-1])
    corr.loc[tissue, '0-1'] = round(spearmanr(parameters[0], parameters[1]).statistic, 3)
    corr.loc[tissue, '1-2'] = round(spearmanr(parameters[1], parameters[2]).statistic, 3)
    corr.loc[tissue, '1-3'] = round(spearmanr(parameters[1], parameters[3]).statistic, 3)
    corr.loc[tissue, '2-3'] = round(spearmanr(parameters[2], parameters[3]).statistic, 3)


corr.to_csv('combined plots/' + centrality + '_subsampled tissues.csv', header = ['s'] + labels)

'''
# %%% Plots

plt.figure()
for c, col in enumerate(corr.columns[1:]):
    plt.scatter(corr['s'], corr[col], label=labels[c], marker=shapes[c+1], color=colors[c+1])
plt.xlim(0,500)
plt.xticks([73, 237])
plt.legend()
for i, tissue in enumerate(corr.index[::-1]):
    xytext = (100+(i//2)*10, 0+(i//2)*0.15) if corr.loc[tissue, 's'] == 73 else (300+(i//2), 0.7+(i//2)*0.15)
    for col in corr.columns:
        plt.annotate(tissue, (corr.loc[tissue, 's'], corr.loc[tissue, col]), xytext=xytext, color=annot_col[i//2], arrowprops=dict(arrowstyle="->", color=annot_col[i//2]))
plt.xlabel('Sample Size ($s$)')
plt.ylabel('Correlation Coefficient ($r$)')
plt.title('Spearman Correlation among Centrality Metrics \n (for subsampled datasets)')
plt.savefig('combined plots/' + centrality + '_subsampled tissues.pdf')
plt.close()
'''