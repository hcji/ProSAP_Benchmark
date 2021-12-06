# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 09:23:50 2021

@author: jihon
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

veh_1 = pd.read_csv('Data/Sav_STS/Vehicle_1.csv')
veh_2 = pd.read_csv('Data/Sav_STS/Vehicle_2.csv')
cas_1 = pd.read_csv('Data/Sav_STS/Stauro_1.csv')
cas_2 = pd.read_csv('Data/Sav_STS/Stauro_2.csv')

accs = cas_1['Accession'].values
cols = cas_1.columns[1:11]
corrs_1, corrs_2 = [], []
for a in accs:
    c1 = cas_1.loc[cas_1['Accession'] == a, cols].values
    c2 = cas_2.loc[cas_2['Accession'] == a, cols].values
    v1 = veh_1.loc[veh_1['Accession'] == a, cols].values
    v2 = veh_2.loc[veh_2['Accession'] == a, cols].values
    
    if min(len(c1), len(c2), len(v1), len(v2)) > 0:
        cor_1 = pearsonr(c1[0,:], c2[0,:])[0]
        cor_2 = pearsonr(v1[0,:], v2[0,:])[0]
        corrs_1.append(cor_1)
        corrs_2.append(cor_2)


plt.figure(figsize = (6, 8), dpi = 300)
plt.subplot(2,1,1)
sns.distplot(corrs_1, bins = 200, hist=True, color = 'r')
plt.xlabel('Pearson correlation', fontsize = 17)
plt.ylabel('Density', fontsize = 17)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlim((0.7, 1))
    
plt.subplot(2,1,2)
sns.distplot(corrs_2, bins = 200, hist=True, color = 'g')
plt.xlabel('Pearson correlation', fontsize = 17)
plt.ylabel('Density', fontsize = 17)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.xlim((0.7, 1))

