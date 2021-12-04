# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 10:15:29 2021

@author: jihon
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
# from adjustText import adjust_text

true_tar = ['Q64133','Q8BW75']
drug_dir = 'Result/Harmine_iTSA'

ttest = pd.read_csv('{}/t_test.csv'.format(drug_dir))
limma = pd.read_csv('{}/limma.csv'.format(drug_dir))
edger = pd.read_csv('{}/edgeR.csv'.format(drug_dir))
deseq = pd.read_csv('{}/DESeq2.csv'.format(drug_dir))

def plot_ranking(tpp, true_tar):
    w = np.array([i for i in tpp.index if tpp.loc[i, 'Accession'] in true_tar])
    plt.scatter(1 + np.arange(len(tpp)), tpp.loc[:,'-logAdjPval'], s = 25)
    plt.scatter(1 + w, tpp.loc[w,'-logAdjPval'], color='r', s = 25)
    '''
    texts = []
    for i, s in enumerate(true_tar):
        x = w[i]
        y = tpp.loc[x,'Score']
        texts.append(plt.text(x+200, y, s, fontsize = 14))
    adjust_text(texts, force_points=0.01, force_text=0.01)
    '''
    plt.xlabel('rank', fontsize = 17)
    plt.xticks(fontsize = 14, rotation=45)
    plt.yticks(fontsize = 14)


plt.figure(figsize = (20, 4), dpi = 300)
plt.subplot(141)
plt.ylabel('-log adjusted p-value', fontsize = 17)
plot_ranking(ttest, true_tar)
plt.subplot(142)
plot_ranking(limma, true_tar)
plt.subplot(143)
plot_ranking(edger, true_tar)
plt.subplot(144)
plot_ranking(deseq, true_tar)




