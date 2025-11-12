# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 20:57:57 2025

@author: HP
"""

# %% Packages
import numpy as np
from corals.correlation.full.matmul import full_matmul_symmetrical as correlation
from corals.correlation.utils import derive_pvalues
from statsmodels.stats.multitest import multipletests as fdr
import pandas as pd
import random
import matplotlib.pyplot as plt
from time import time

# %% Seed
import os
seed_file = 'seed.txt' 
if not os.path.exists(seed_file):
    seed = int(time())
    with open(seed_file, 'w') as file:
        print(seed, file = file)
    file.close()
else:
    with open(seed_file, 'r') as file:
        seed = int(file.read())
np.random.seed(seed)
random.seed(seed)


# %% Co-expression Network Construction
def compute_deg(data):
    r = correlation(data, correlation_type='spearman')
    p = derive_pvalues(r, len(data))
    
    indices = np.triu_indices(n, k= 1)
    p_flat = p[indices]
    rejected, _, _, _ = fdr(p_flat, alpha=0.01, method='fdr_bh')
    
    adj_matrix = np.zeros((n,n))
    adj_matrix[indices] = rejected
    adj_matrix = adj_matrix + adj_matrix.T   
    return np.count_nonzero(adj_matrix[0, :])


# %% Population Dataset Generation
print('Population Dataset Generation')
n = 201
B = 1000
pop_size = 1000000

data_pop = pd.DataFrame(index=range(pop_size), columns=range(n))
for i in range(1,n):
    data_pop[i] = np.random.normal(size = pop_size)

# %% For all s -- the main funciton
for s in [670, 237, 73]:
    # %%% All Computations with 1% FDR (Tasks Estimator Bias 7, Bootstrap Bias 1 to 4)
    print("Estimator and Bootstrap Bias Estimation \t ", s)

    # Biases
    bias_estimator = pd.DataFrame(columns=['deg_pop', 'bias'])
    bias_bootstrap = pd.DataFrame(columns=['deg_obs', 'bias_p1', 'bias_p2', 'bias_p3'])
    bias_test = pd.DataFrame(columns=['deg_pop', 'bias'])

    for deg in range(1, n):
        print(deg)
        
        print('Population Data')
        chosen_vars = random.sample(range(1,n), k=deg)
        alpha = np.zeros(shape=(n,))
        alpha[chosen_vars] = 0.995
        alpha[0] = 0.1
        beta = 0.01
        
        data_pop.loc[:,0] = [1]*pop_size
        noise = np.random.normal(size = pop_size)
        
        x = []
        for i in range(0, pop_size, 10000):
            x.append(data_pop.loc[i:(i+10000-1), :] * alpha)
            
        y = []
        for i in range(len(x)):
            y.append(x[i].sum(axis = 1))
        
        z = pd.concat(y, ignore_index=False)
        data_pop[0] = z + (beta * noise)
        deg_pop = compute_deg(data_pop)
        
        # Data Sample
        print('Observed Data Sampling')
        deg_b = np.zeros(shape = (B,), dtype = int)
        for b in range(B):
            data_sample = data_pop.loc[random.sample(range(pop_size), k=s), :]
            deg_b[b] = compute_deg(data_sample)
        deg_obs = np.mean(deg_b)
        
        bias_estimator.loc[deg,:] = [deg_pop, deg_obs-deg_pop]
        bias_bootstrap.loc[deg, 'deg_obs'] = deg_b[-1]
        
        # Bootstrap Sample
        print("Bootstrapping")
        deg_b = np.zeros(shape = (B,), dtype = int)
        for b in range(B):
            bootstrap_sample = data_sample.loc[random.choices(data_sample.index, k = s), :]
            deg_b[b] = compute_deg(bootstrap_sample)
        deg_bootstrap = np.mean(deg_b)
        std_bootstrap = np.std(deg_b)
        
        bias_bootstrap.loc[deg, 'bias_1'] = deg_bootstrap - bias_bootstrap.loc[deg, 'deg_obs']
        bias_bootstrap.loc[deg, 'bias_2'] = bias_bootstrap.loc[deg, 'bias_1'] - (2 * std_bootstrap)
        bias_bootstrap.loc[deg, 'bias_3'] = bias_bootstrap.loc[deg, 'bias_1'] - std_bootstrap
        
        bias_test.loc[deg,:] = [deg_pop, deg_bootstrap - deg_pop]
    
    # %%% PLOTS
    print('Plotting')
    plt.figure()
    plt.scatter(bias_estimator['deg_pop'], bias_estimator['bias'])
    plt.axline((0,0), slope=0, color='#000000')
    plt.xlabel('Population Degree')
    plt.ylabel('Bias \n E[Observed Degree] - Actual Degree')
    plt.title('Estimator Bias -- 1% FDR')
    plt.savefig('s_' + str(s) + '/Estimator Bias FDR 0.01.pdf')
    plt.close()

    plt.figure()
    plt.scatter(bias_bootstrap['deg_obs'], bias_bootstrap['bias_1'])
    plt.axline((0,0), slope=0, color='#000000')
    plt.xlabel('Observed Degree of the 1000th sample')
    plt.ylabel('Bias \n Bootstrap Degree - Observed Degree')
    plt.title('Bootstrap Bias')
    plt.savefig('s_' + str(s) + '/Bootstrap Bias p1.pdf')
    plt.close()
