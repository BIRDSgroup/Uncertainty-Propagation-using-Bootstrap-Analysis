# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 21:53:09 2024

@author: HP
"""
# %% Tasks
'''
Estimator Bias
1. r cutoff 0.1
2. r cutoff 0.5
3. r cutoff 0.8
4. p cutoff 0.01
5. p cutoff 0.05
6. 10% FDR cutoff (0.1)

7. 1% FDR cutoff (0.01) 
1. Bootstrap Bias mu - deg_obs v/s deg_obs
2. Bootstrap Bias w.r.t Actual Degree
3. mu - 2sigma - deg_obs v/s deg_obs
4. mu - sigma - deg_obs v/s deg_obs

'''

# %% Packages
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests as fdr
import pandas as pd
import random
import matplotlib.pyplot as plt

# %% Co-expression Network Construction
def compute_deg(data):
    p = spearmanr(data)[1]
    indices = np.triu_indices(n, k= 1)
    p_flat = p[indices]
    rejected, _, _, _ = fdr(p_flat, alpha=0.01, method='fdr_bh')
    
    adj_matrix = np.zeros((n,n))
    adj_matrix[indices] = rejected
    adj_matrix = adj_matrix + adj_matrix.T   
    return np.count_nonzero(adj_matrix[0, :])

# %% PART -- 1 = All Computations with 1% FDR (Tasks Estimator Bias 7, Bootstrap Bias 1 to 4)
n = 100
s = 73
B = 1000

# Biases
bias_estimator = pd.DataFrame(columns=['deg_pop', 'bias'])
bias_bootstrap = pd.DataFrame(columns=['deg_obs', 'bias_p1', 'bias_p2', 'bias_p3'])
bias_test = pd.DataFrame(columns=['deg_pop', 'bias'])


for deg in range(1, n):
    print(deg)
    
    # Population
    print('Population Dataset Generation')
    chosen_vars = random.sample(range(1,n), k=deg)
    alpha = np.zeros(shape=(n,))
    alpha[chosen_vars] = 0.995
    alpha[0] = 0.1
    beta = 0.01
    pop_size = 1000000
    
    data_pop = pd.DataFrame(index=range(pop_size), columns=range(n))
    for i in range(1,n):
        data_pop[i] = np.random.normal(size = pop_size)
    data_pop[0] = pd.Series([1]*pop_size)
    noise = np.random.normal(size = pop_size)
    data_pop[0] = (data_pop @ alpha) + (beta * noise)
    deg_pop = compute_deg(data_pop)
    
    # Data Sample
    print('Observed Data Sampling')
    deg_b = np.zeros(shape = (B,), dtype = int)
    for b in range(B):
        data_sample = data_pop.loc[random.sample(list(data_pop.index), k=s), :]
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
    
# %% PLOTS
plt.figure()
plt.scatter(bias_estimator['deg_pop'], bias_estimator['bias'])
plt.xlabel('Population / Actual Degree')
plt.ylabel('Bias \n E[Observed Degree] - Actual Degree')
plt.title('Estimator Bias -- 1% FDR')
plt.savefig('s ' + str(s) + '/Part 1/Estimator Bias 1% FDR.png')

plt.figure()
plt.scatter(bias_bootstrap['deg_obs'], bias_bootstrap['bias_1'])
plt.xlabel('Observed Degree of the 1000th sample')
plt.ylabel('Bias \n Bootstrap Degree - Observed Degree')
plt.title('Bootstrap Bias')
plt.savefig('s ' + str(s) + '/Part 1/Bootstrap Bias -- p1.png')

plt.figure()
plt.scatter(bias_test['deg_pop'], bias_test['bias'])
plt.xlabel('Population / Actual Degree')
plt.ylabel('Bias \n Bootstrap Degree - Actual Degree')
plt.title('Bootstrap Bias w.r.t Actual Degree')
plt.savefig('s ' + str(s) + '/Part 1/Bootstrap Bias #2.png')

plt.figure()
plt.scatter(bias_bootstrap['deg_obs'], bias_bootstrap['bias_2'])
plt.xlabel('Observed Degree of the 1000th sample')
plt.ylabel('Bias \n (Bootstrap Degree - 2 * Sigma) - Observed Degree')
plt.title('Bias of mu - 2sigma')
plt.savefig('s ' + str(s) + '/Part 1/Bootstrap Bias -- p2.png')

plt.figure()
plt.scatter(bias_bootstrap['deg_obs'], bias_bootstrap['bias_3'])
plt.xlabel('Observed Degree of the 1000th sample')
plt.ylabel('Bias \n (Bootstrap Degree - Sigma) - Observed Degree')
plt.title('Bias of mu - sigma')
plt.savefig('s ' + str(s) + '/Part 1/Bootstrap Bias -- p3.png')

# %%
%reset -f # this code is for spyder only -- this is not a python command but a spyder / jupyter command to delete all local variable.
# locals().clear()
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests as fdr
import pandas as pd
import random
import matplotlib.pyplot as plt


# %% PART -- 2
def compute_deg_test(data):
    r, p = spearmanr(data)
    indices = np.triu_indices(n, k= 1)
    p_flat = p[indices]
    rejected, _, _, _ = fdr(p_flat, alpha=0.1, method='fdr_bh')
    
    adj_matrix = np.zeros((n,n))
    adj_matrix[indices] = rejected
    adj_matrix = adj_matrix + adj_matrix.T
    degs = (np.count_nonzero(abs(r[0,:]) >= 0.1), np.count_nonzero(abs(r[0,:]) >= 0.5), np.count_nonzero(abs(r[0,:]) > 0.8), np.count_nonzero(p[0,:] <= 0.01), np.count_nonzero(p[0,:] <= 0.05), np.count_nonzero(adj_matrix[0,:]))
    return degs
    

n = 100
s = 73
B = 1000

# Biases
col_names = [['deg_pop_' + str(i), 'bias_' + str(i)] for i in range(1,7)]
col_names = sum(col_names, [])
bias_estimator = pd.DataFrame(columns= col_names)

for deg in range(1, n):
    print(deg)
    
    # Population
    print('Population Dataset Generation')
    chosen_vars = random.sample(range(1,n), k=deg)
    alpha = np.zeros(shape=(n,))
    alpha[chosen_vars] = 0.995
    alpha[0] = 0.1
    beta = 0.01
    pop_size = 1000000
    
    data_pop = pd.DataFrame(index=range(pop_size), columns=range(n))
    for i in range(1,n):
        data_pop[i] = np.random.normal(size = pop_size)
    data_pop[0] = pd.Series([1]*pop_size)
    noise = np.random.normal(size = pop_size)
    data_pop[0] = (data_pop @ alpha) + (beta * noise)
    deg_pop = compute_deg_test(data_pop)
    
    # Data Sample
    print('Observed Data Sampling')
    deg_b = np.zeros(shape = (B,6), dtype = int)
    for b in range(B):
        data_sample = data_pop.loc[random.sample(list(data_pop.index), k=s), :]
        deg_b[b,:] = compute_deg_test(data_sample)
    deg_obs = np.mean(deg_b, axis=0)
    
    index = len(bias_estimator)
    for i in range(1,7):
        bias_estimator.loc[index, ['deg_pop_' + str(i), 'bias_' + str(i)]] = [deg_pop[i-1], deg_obs[i-1] - deg_pop[i-1]]
    


titles = ['Estimator Bias -- r value cutoff 0.1', 'Estimator Bias -- r value cutoff 0.5', 'Estimator Bias -- r value cutoff 0.8', 'Estimator Bias -- p value cutoff 0.01', 'Estimator Bias -- p value cutoff 0.05', 'Estimator Bias -- 10% FDR']
filename = ['Estimator Bias r 0.1.png', 'Estimator Bias r 0.5.png', 'Estimator Bias r 0.8.png', 'Estimator Bias p 0.01.png', 'Estimator Bias p 0.05.png', 'Estimator Bias 10% FDR.png']
for i in range(1, 7):
    plt.figure()
    plt.scatter(bias_estimator['deg_pop_' + str(i)], bias_estimator['bias_' + str(i)])
    plt.xlabel('Population / Actual Degree')
    plt.ylabel('Bias \n E[Observed Degree] - Actual Degree')
    plt.title(titles[i-1])
    plt.savefig('s ' + str(s) + '/Part 2/' +  filename[i-1])
    