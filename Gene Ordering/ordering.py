# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 22:44:55 2025

@author: HP
"""

import numpy as np
import pandas as pd
import random
import os
from time import time

#Tissues = ['Kidney_Cortex', 'Lung', 'Muscle_Skeletal', 'Pancreas', 'Pituitary', 'Stomach', 'Thyroid', 'Whole_Blood', 'Vagina']

Tissues=["Whole_Blood", "Muscle_Skeletal", "Lung", "Skin_Sun_Exposed_Lower_leg", "Thyroid", "Pancreas", "Brain_Cortex", "Pituitary", "Brain_Cerebellum", "Stomach", "Brain_Caudate_basal_ganglia", "Kidney_Cortex", "Brain_Substantia_nigra", "Uterus", "Vagina", "Brain_Amygdala"]

# %% Set Seed
seed_file = 'seed.txt'
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

# %% Open Gencode
gencode = pd.read_csv('gencode.csv')
order = random.sample(range(len(gencode)), len(gencode))
ordered_gencode = gencode.loc[order,:]

# %% Tissue List
for tissue in Tissues:
    genes = pd.read_csv('../Original Dataset/Preprocessed Files/' + tissue + '/genes.csv')
    index_map = pd.Series(genes.index, index=genes['gene_id'])
    
    df = ordered_gencode[ordered_gencode['gene_id'].isin(genes['gene_id'])]
    tissue_order = index_map[df['gene_id']]
    tissue_order.to_csv(tissue + '_order.csv', index = False, header=False)