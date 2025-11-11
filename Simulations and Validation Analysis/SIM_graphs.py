# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 14:45:55 2024

@author: HP
"""

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests as fdr
import random

from time import time

import matplotlib.pyplot as plt
import seaborn as sns

import plots

# %% Global Variables / Inputs
n = 400
B = 1000
klim = 100

SampleSizes = [73, 237, 670]
'''
For reference = 706, SampleSizes = [53, 73, 237, 670]
For reference = 1500, SampleSizes = [53, 73, 237, 670, 706, 1000]
'''
s = 706        # Sample Size of the Reference Dataset


# %% Co-expression Network Construction
def compute_deg(data):
    p = spearmanr(data)[1]
    indices = np.triu_indices(n, k=1)
    p_flat = p[indices]
    rejected, _, _, _ = fdr(p_flat, alpha=0.01, method='fdr_bh')
    
    adj_matrix = np.zeros((n,n))
    adj_matrix[indices] = rejected
    adj_matrix = adj_matrix + adj_matrix.T
        
    plt.figure()
    ax = sns.heatmap(adj_matrix[n//2:, :n//2], cmap='Reds', vmin=0, vmax=1, xticklabels=False, yticklabels=False)
    plt.title('Adjacency Matrix \n s = ' + str(len(data)))
    plt.xlabel('Genes - Set $S_1$', fontsize=13)
    plt.ylabel('Genes - Set $S_2$', fontsize=13)
    plt.tight_layout()
    plt.savefig('Simulated Data/graphs/heatmap_' + str(len(data)) + '.pdf')
    plt.close()
    
    deg = np.sum(adj_matrix[n//2:, :n//2], axis = 1) # 3rd Quadrant of the adjacency matrix
    
    plt.figure()
    plt.hist(deg)
    plt.title('Degree Distribution \n s = ' + str(len(data)))
    plt.xlabel('Degrees', fontsize=13)
    plt.ylabel('Frequency', fontsize=13)
    plt.tight_layout()
    plt.savefig('Simulated Data/graphs/degreeDist_' + str(len(data)) + '.pdf')
    plt.close()
    
    return deg

def Bootstrapping(data_sample):
    # data_sample.shape() = s x n
    print("Bootstrapping")
    deg_b = np.zeros(shape = (B,n//2), dtype = int)
    for b in range(B):
        if b % 100 == 0:
            print(b)
        bootstrap_sample = data_sample.loc[random.choices(data_sample.index, k = s),:]
        #deg_b[b,:] = compute_deg(bootstrap_sample)
    deg_bootstrap = np.mean(deg_b, axis = 0)
    std_bootstrap = np.std(deg_b, axis = 0)
    return deg_bootstrap, std_bootstrap


# %% Seed
import os
seed_file = 'Simulated Data/seed.txt'
if not os.path.exists(seed_file):
    seed = int(time())
    with open(seed_file, 'w') as file:
        print(seed, file = file)
    file.close()
else:
    with open(seed_file, 'r') as file:
        seed = int(file.read())
random.seed(seed)
np.random.seed(seed)

# Random Ordering
order = random.sample(range(n//2), n//2)

# %% Population Dataset Generation
print('Population Dataset Generation')
pop_size = 1000000

GT = []
data_pop = pd.DataFrame(index=range(pop_size), columns=range(n))
for i in range(n//2):
    data_pop[i] = np.random.normal(size = pop_size)

for i in range(n//2, n):
    deg = random.randint(1, n//2)
    chosen_vars = random.sample(range(n//2), k=deg)
    alpha = np.zeros(shape=(n//2 + 1,))
    alpha[chosen_vars] = 0.995
    alpha[-1] = 0.1
    beta = 0.01
    
    GT.append(deg)
    
    data_pop[i] = pd.Series([1]*pop_size)
    noise = np.random.normal(size = pop_size)
    #data_pop[i] = data_pop.loc[:,range(n//2)] @ alpha[:-1] + (data_pop.loc[:,i] * alpha[-1]) + (beta * noise)
    
    data_pop.loc[:,i] = (data_pop.loc[:,i] * alpha[-1]) + (beta * noise)
    for j in range(n//2):
        data_pop.loc[:,i] += (data_pop.loc[:,j] * alpha[j])
    
dataset = data_pop.copy() ## dataset.shape() = s x n
deg_pop_orig = compute_deg(dataset)
del data_pop

plots.DegComp(deg_pop_orig, GT, s=None, pop=True)
print('Avg Degree of Population = ', np.mean(deg_pop_orig))

# %% Observed and Bootstrap Dataset -- I (used for validation and replication analysis only.)
# Data Sample
print('Observed Data Sampling')
data_sample = dataset.loc[random.sample(list(dataset.index), k=s), :]
deg_obs_1 = compute_deg(data_sample)

# Bootstrapping
deg_bootstrap_1, std_bootstrap_1 = Bootstrapping(data_sample)

plots.DegComp(deg_obs_1, deg_pop_orig, s=s, pop=False)
print('Avg Degree of', s, '=', np.mean(deg_obs_1), sep = ' ')

# %% Dataset -- II
for s in SampleSizes:   
    # %%% Observed and Bootstrap Dataset -- II    
    # Data Sample
    print('Observed Data Sampling -- ', s, sep = '')
    data_sample = dataset.loc[random.sample(list(dataset.index), k=s), :]
    deg_obs_2 = compute_deg(data_sample)
    
    # Bootstrapping
    deg_bootstrap_2, std_bootstrap_2 = Bootstrapping(data_sample)
    
    plots.DegComp(deg_obs_2, deg_pop_orig, s=s, pop=False)
    print('Avg Degree of', s, '=', np.mean(deg_obs_2), sep = ' ')

