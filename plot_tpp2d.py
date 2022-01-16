# -*- coding: utf-8 -*-
"""
Created on Sun Dec  5 07:46:49 2021

@author: jihon
"""


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


data = pd.read_csv('Data/TPP2D/panobinostat_tpp2d_cell.csv')


def plot_tpp2d(data, ProteinAccession):
    pltdata = data[data['clustername'] == ProteinAccession]
    conc = np.unique(pltdata['conc'])
    temp = np.unique(pltdata['temperature'])
    img = np.zeros((len(conc), len(temp)))
    for i in pltdata.index:
        a = np.where(conc == pltdata.loc[i,'conc'])[0][0]
        b = np.where(temp == pltdata.loc[i,'temperature'])[0][0]
        img[a, b] = pltdata.loc[i,'rel_value']
    img = pd.DataFrame(img)
    img.index = conc
    img.columns = temp
    
    sns.heatmap(img, cbar=True)
    plt.xlabel('temperture', fontsize = 14)
    plt.ylabel('drug concentration', fontsize = 14)


plt.figure(figsize = (15, 4), dpi = 300)
plt.subplot(131)
plot_tpp2d(data, 'HDAC1')
plt.subplot(132)
plot_tpp2d(data, 'HDAC2')
plt.subplot(133)
plot_tpp2d(data, 'TTC38')

