# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 14:06:22 2025

@author: HP

Working Directory is: Replication Analysis 2\Original Dataset\Muscle_Skeletal_recount
"""

import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

cov_data = pd.read_csv('covariates.csv', header = 0, index_col=0)
df = pd.read_csv('Muscle_Skeletal_recount normalized.csv', header=0, index_col=0)
corr_1 = pd.Series(index = df.index)
for idx in df.index:
    print(idx)
    corr_1[idx] = spearmanr(df.loc[idx,:], cov_data['Age']).statistic
    
df = pd.read_csv('Muscle_Skeletal_recount.csv', header=0, index_col=0)
corr_2 = pd.Series(index = df.index)
for idx in df.index:
    print(idx)
    corr_2[idx] = spearmanr(df.loc[idx,:], cov_data['Age']).statistic

plt.figure()
plt.scatter(corr_1, corr_2, color='#CC79A7')
plt.axline((0,0), slope=1, color='#000000')
plt.xlabel('Without Covariate Adjustment')
plt.ylabel('With Covariate Adjustment')

print(corr_1.nlargest(n=20).index)
print(corr_2.nlargest(n=20).index)

