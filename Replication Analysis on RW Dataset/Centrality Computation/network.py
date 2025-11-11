#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 22:02:38 2024

@author: sugyani
"""

from envVar import setThreads
setThreads(1)

# %% Packages
import pandas as pd
from numpy import round
import sys

import bootstrap
import centrality

# %% Input
arg = sys.argv
tissue = arg[1]

folder = tissue
if '_' in tissue and tissue.split(sep='_')[-1].isdigit():
    folder = "_".join(tissue.split(sep='_')[:-1])
in_file = '../Original Dataset/Preprocessed Files/' + folder + '/' + tissue + '.csv'

# %% 
gexp = pd.read_csv(in_file, sep = ',', header=0, index_col=0)
n, s = gexp.shape

edges = bootstrap.CoExpressionNetwork(gexp.T)

deg = centrality.DegreeCentrality(edges, n)
with open(tissue + '/original_degree.csv', 'w') as file:
    print(*deg, sep = ',', end = '\n', file=file)
file.close()

pr = centrality.PageRankCentrality(edges, n)
with open(tissue + '/original_pagerank.csv', 'w') as file:
    print(*round(pr,3), sep=',', end='\n', file=file)
file.close()
