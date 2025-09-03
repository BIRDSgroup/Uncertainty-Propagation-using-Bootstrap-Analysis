#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 21:13:36 2024

@author: sugyani
"""

from envVar import setThreads
setThreads(1)

# %% Packages
import pandas as pd
from numpy import round
import random
import sys

import bootstrap
import centrality

# %% Input
arg = sys.argv
tissue = arg[1]
b = int(arg[2])

folder = tissue
if '_' in tissue and tissue.split(sep='_')[-1].isdigit():
    folder = "_".join(tissue.split(sep='_')[:-1])
in_file = '../Original Dataset/Preprocessed Files/' + folder + '/' + tissue + '.csv'

deg_file = 'degree/sample_degree_' + str(b) + '.csv'
pr_file = 'pagerank/sample_pagerank_' + str(b) + '.csv'

# %% Set the Seeds
with open('seeds/seed_' + str(b) + '', 'r') as seed:
    random.seed(seed.read())
seed.close()

# %% Bootstrapping
gexp = pd.read_csv(in_file, sep = ',', header=0, index_col=0)
n, s = gexp.shape
edges = bootstrap.BootstrapSample(gexp.T)

# %% Degree Centrality
deg = centrality.DegreeCentrality(edges, n)
with open(deg_file, 'w') as file:
    print(*deg, sep=',', end='\n', file=file)
file.close()

# %% Pagerank Centrality
pr = centrality.PageRankCentrality(edges, n)
with open(pr_file, 'w') as file:
    print(*round(pr,3), sep=',', end='\n', file=file)
file.close()

