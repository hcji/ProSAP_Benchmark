# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 15:24:50 2021

@author: jihon
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('Data/Sav_STS/Staturosporine_TPP_data_Savitski.csv')
kins = data.loc[data.loc[:,'is kinase?']=='YES','Accession'].values

tpp = pd.read_csv('Result/Sav_STS/TPP.csv')
tpp_nofilter = pd.read_csv('Result/Sav_STS/TPP_no_filter.csv')

wh = np.where([p in list(kins) for p in tpp_nofilter['Accession']])[0]

plt.figure(figsize=(6, 4.5), dpi=300)
plt.plot(tpp_nofilter['Rep1delta_Tm'], tpp_nofilter['Rep2delta_Tm'], '.' , label='Others')
plt.plot(tpp_nofilter.loc[wh, 'Rep1delta_Tm'], tpp_nofilter.loc[wh, 'Rep2delta_Tm'], '.', label='Kinase', color='#C00000')
plt.plot([-10, 18], [-10, 18], color = 'black', linestyle='--')
plt.plot([0, 0], [-10, 18], color='gray', linestyle='--')
plt.plot([-10, 18], [0, 0], color='gray', linestyle='--')
plt.xlabel('Delta Tm of Replicate 1')
plt.ylabel('Delta Tm of Replicate 2')
plt.legend()

tpp_rep1 = tpp_nofilter.sort_values(by='Rep1pVal (-log10)', ascending=False)
tpp_rep1 = tpp_rep1.reset_index()
tpp_rep2 = tpp_nofilter.sort_values(by='Rep2pVal (-log10)', ascending=False)
tpp_rep2 = tpp_rep2.reset_index()

tpp_ifkin = np.cumsum([i in kins for i in tpp['Accession']])
tpp_nofilter_ifkin = np.cumsum([i in kins for i in tpp_nofilter['Accession']])
tpp_rep1_ifkin = np.cumsum([i in kins for i in tpp_rep1['Accession']])
tpp_rep2_ifkin = np.cumsum([i in kins for i in tpp_rep2['Accession']])

plt.figure(figsize=(6,4.5), dpi=300)
plt.plot(np.arange(1, 101) - tpp_ifkin[:100], tpp_ifkin[:100], label='TPP with filter')
plt.plot(np.arange(1, 101) - tpp_nofilter_ifkin[:100], tpp_nofilter_ifkin[:100], label='TPP without filter')
plt.plot(np.arange(1, 101) - tpp_rep1_ifkin[:100], tpp_rep1_ifkin[:100], label='TPP replicate 1')
plt.plot(np.arange(1, 101) - tpp_rep2_ifkin[:100], tpp_rep2_ifkin[:100], label='TPP replicate 2')
plt.xlabel('Number of false positive')
plt.ylabel('Number of true positive')
plt.legend()

