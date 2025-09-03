# -*- coding: utf-8 -*-
"""
Created on Mon Sep  1 15:27:55 2025

@author: HP
"""

import pandas as pd
metadata = pd.read_csv('covariates.txt', sep='\t', header=None, index_col=0).T
cov_df = pd.DataFrame({
    'Age': metadata.iloc[:,10],
    'Gender': metadata.iloc[:,11]
    })
cov_df.index=metadata.iloc[:,0] + '_COUNT'
cov_df['Age'] = cov_df['Age'].str.split().str[-1]
cov_df['Gender'] = cov_df['Gender'].str.split().str[-1]
cov_df['Gender'] = (cov_df['Gender'] == 'M').astype(int)
df = pd.read_csv('Muscle_Skeletal_recount normalized.csv', header=0, index_col=0)
cov_df = cov_df.loc[df.columns,:]

cov_df.to_csv('covariates.csv', header = True, index=True)