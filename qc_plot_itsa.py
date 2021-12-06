# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 10:38:11 2021

@author: jihon
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

data_1 = pd.read_excel('Data/MTX_PISA/MTX_PISA_lysate.xlsx')
data_2 = pd.read_excel('Data/5_FU_PISA/5_FU_PISA.xlsx')
data_3 = pd.read_excel('Data/Harmine_iTSA/Harmine.xlsx')
data_3.iloc[:,1:] = 2 ** data_3.iloc[:,1:]
data_4 = pd.read_csv('Data/Ball_STS_iTSA/Staturosporine_iTSA_data_Ball.csv')
data_4.iloc[:,1:] = 2 ** data_4.iloc[:,1:]

def cal_rsd(data_1):
    stds_1, stds_2 = [], []
    for i in range(len(data_1)):
        y1 = data_1.iloc[i, 1:6].values
        y2 = data_1.iloc[i, 6:].values
        stds_1.append(np.std(y1) / np.mean(y1))
        stds_2.append(np.std(y2) / np.mean(y2))

    plt.figure(figsize = (6, 8), dpi = 300)
    plt.subplot(2,1,1)
    sns.distplot(stds_1, bins = 200, hist=True, color = 'r')
    plt.xlabel('RSD', fontsize = 17)
    plt.ylabel('Density', fontsize = 17)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlim((0, 0.5))        

    plt.subplot(2,1,2)
    sns.distplot(stds_2, bins = 200, hist=True, color = 'g')
    plt.xlabel('RSD', fontsize = 17)
    plt.ylabel('Density', fontsize = 17)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.xlim((0, 0.5))   

cal_rsd(data_1)
cal_rsd(data_2)
cal_rsd(data_3)
cal_rsd(data_4)

