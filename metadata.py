# -*- coding: utf-8 -*-
"""
Created on Sat Mar 29 09:52:17 2025

@author: HP
"""


import pandas as pd

Tissues = ['Brain_Amygdala', 'Brain_Caudate_basal_ganglia', 'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Substantia_nigra', 'Kidney_Cortex', 'Lung', 'Muscle_Skeletal', 'Pancreas', 'Pituitary', 'Skin_Sun_Exposed_Lower_leg', 'Stomach', 'Thyroid', 'Uterus', 'Vagina', 'Whole_Blood', 'Muscle_Skeletal_73', 'Muscle_Skeletal_237', 'Whole_Blood_73', 'Whole_Blood_237', 'Lung_237', 'Lung_73']

for tissue in Tissues:
    folder = tissue
    if '_' in tissue and tissue.split(sep='_')[-1].isdigit():
        folder = "_".join(tissue.split(sep='_')[:-1])
    in_file = '../Original Dataset/Preprocessed Files/' + folder + '/' + tissue + '.csv'
    
    gexp = pd.read_csv(in_file, sep = ',', header=0, index_col=0)
    n, s = gexp.shape
    print(tissue, n, s, sep = '\t')
    
