# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 14:22:22 2021

@author: jihon
"""


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text

true_tar = ['TTC38', 'HDAC1', 'HDAC2']

tpp2d = pd.read_csv('Result/TPP2D/panobinostat_tpp2d_cell.csv')

def plot_ranking(tpp, true_tar):
    w = np.array([i for i in tpp.index if tpp.loc[i, 'clustername'] in true_tar])
    plt.scatter(1 + np.arange(len(tpp)), tpp.loc[:,'F_statistic'], s = 25)
    plt.scatter(1 + w, tpp.loc[w,'F_statistic'], color='r', s = 25)
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
    
plt.figure(figsize = (8, 6), dpi = 300)
plot_ranking(tpp2d, true_tar)

