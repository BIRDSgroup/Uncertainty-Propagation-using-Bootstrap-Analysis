#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 15:39:13 2023

@author: sugyani
"""

# %% Packages
import pandas as pd
import sys

# %% Input
arg = sys.argv
tissue = arg[1]
B = int(arg[2])

# %% Merge Degree
with open(tissue + '/sample_degree.csv', 'w') as file:
    for b in range(B):
        deg = pd.read_csv('degree/sample_degree_' + str(b) + '.csv', header=None, index_col=None, sep=',').iloc[0,:]
        print(b, *deg, sep=',', end='\n', file=file)
file.close()


# %% Merge Pagerank
with open(tissue + '/sample_pagerank.csv', 'w') as file:
    for b in range(B):
        pr = pd.read_csv('pagerank/sample_pagerank_' + str(b) + '.csv', header=None, index_col=None, sep=',').iloc[0,:]
        print(b, *pr, sep=',', end='\n', file=file)
file.close()