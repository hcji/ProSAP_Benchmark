# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 10:15:29 2021

@author: jihon
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text

true_tar = ['TTC38', 'HDAC6', 'HDAC8']
drug_dir = 'Result/Panobinostat'

tpp = pd.read_csv('{}/tpp.csv'.format(drug_dir))
npa = pd.read_csv('{}/NPARC.csv'.format(drug_dir))
euc = pd.read_csv('{}/Euclidean.csv'.format(drug_dir))
che = pd.read_csv('{}/Chebyshev.csv'.format(drug_dir))

def plot_ranking(tpp, true_tar):
    w = np.array([i for i in tpp.index if tpp.loc[i, 'Accession'] in true_tar])
    plt.scatter(1 + np.arange(len(tpp)), tpp.loc[:,'Score'], s = 25)
    plt.scatter(1 + w, tpp.loc[w,'Score'], color='r', s = 25)
    '''
    texts = []
    for i, s in enumerate(true_tar):
        x = w[i]
        y = tpp.loc[x,'Score']
        texts.append(plt.text(x+200, y, s, fontsize = 14))
    adjust_text(texts, force_points=0.01, force_text=0.01)
    '''
    plt.xlabel('rank', fontsize = 17)
    plt.ylabel('score', fontsize = 17)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)


plt.figure(figsize = (20, 4), dpi = 300)
plt.subplot(141)
plot_ranking(tpp, true_tar)
plt.subplot(142)
plot_ranking(npa, true_tar)
plt.subplot(143)
plot_ranking(euc, true_tar)
plt.subplot(144)
plot_ranking(che, true_tar)




