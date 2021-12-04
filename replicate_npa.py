# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 10:21:42 2021

@author: jihon
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('Data/Sav_STS/Staturosporine_TPP_data_Savitski.csv')
kins = data.loc[data.loc[:,'is kinase?']=='YES','Accession'].values

euclidean = pd.read_csv('Result/Sav_STS/Euclidean.csv')

tpp_rep1 = euclidean.sort_values(by='Rep1delta_SH', ascending=False)
tpp_rep1 = tpp_rep1.reset_index()
tpp_rep2 = euclidean.sort_values(by='Rep2delta_SH', ascending=False)
tpp_rep2 = tpp_rep2.reset_index()

euclidean_ifkin = np.cumsum([i in kins for i in euclidean['Accession']])
euclidean_rep1_ifkin = np.cumsum([i in kins for i in tpp_rep1['Accession']])
euclidean_rep2_ifkin = np.cumsum([i in kins for i in tpp_rep2['Accession']])

plt.figure(figsize=(6,4.5), dpi=300)
plt.plot(np.arange(1, 101) - euclidean_ifkin[:100], euclidean_ifkin[:100], label='Euclidean with replicates')
plt.plot(np.arange(1, 101) - euclidean_rep1_ifkin[:100], euclidean_rep1_ifkin[:100], label='Euclidean replicate 1')
plt.plot(np.arange(1, 101) - euclidean_rep2_ifkin[:100], euclidean_rep2_ifkin[:100], label='Euclidean replicate 2')
plt.xlabel('Number of false positive')
plt.ylabel('Number of true positive')
plt.legend()