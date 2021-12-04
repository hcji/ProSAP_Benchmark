# -*- coding: utf-8 -*-
"""
Created on Sat Dec  4 08:53:03 2021

@author: jihon
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
# from adjustText import adjust_text


veh_1 = pd.read_csv('Data/Luminespib/Vehicle_1.csv')
veh_2 = pd.read_csv('Data/Luminespib/Vehicle_2.csv')
cas_1 = pd.read_csv('Data/Luminespib/Luminespib_1.csv')
cas_2 = pd.read_csv('Data/Luminespib/Luminespib_1.csv')


def meltCurve(T, a, b, pl):
    A = 1 - pl
    B = 1 + np.exp(b - a/T)
    return A / B + pl



def plot_curve(veh_1, cas_1, veh_2, cas_2, true_tar):
    c = veh_1.columns[1:10]
    w1 = np.where(veh_1['Accession'].values == true_tar)[0][0]
    w2 = np.where(veh_2['Accession'].values == true_tar)[0][0]
    u1 = np.where(cas_1['Accession'].values == true_tar)[0][0]
    u2 = np.where(cas_2['Accession'].values == true_tar)[0][0]

    x = np.array([float(t.replace('T', '')) for t in c])
    x1 = np.arange(x[0], x[-1], 0.01)
    
    y1 = veh_1.loc[w1, c].values
    y2 = veh_2.loc[w2, c].values
    y3 = cas_1.loc[u1, c].values
    y4 = cas_2.loc[u2, c].values
    
    paras1 = curve_fit(meltCurve, x, y1, bounds=(0, [15000, 250, 0.3]))[0]
    paras2 = curve_fit(meltCurve, x, y2, bounds=(0, [15000, 250, 0.3]))[0]
    paras3 = curve_fit(meltCurve, x, y3, bounds=(0, [15000, 250, 0.3]))[0]
    paras4 = curve_fit(meltCurve, x, y4, bounds=(0, [15000, 250, 0.3]))[0]
        
    yy1 = meltCurve(x1, paras1[0], paras1[1], paras1[2])
    yy2 = meltCurve(x1, paras2[0], paras2[1], paras2[2])
    yy3 = meltCurve(x1, paras3[0], paras3[1], paras3[2])
    yy4 = meltCurve(x1, paras4[0], paras4[1], paras4[2])
    
    plt.figure(figsize = (10, 4), dpi = 300)
    plt.subplot(121)
    plt.scatter(x, y1, label = None)
    plt.scatter(x, y3, label = None)
    plt.plot(x1, yy1, label = 'Vehicle_1')
    plt.plot(x1, yy3, label = 'Drug_1')
    plt.xlabel('temperture', fontsize = 17)
    plt.ylabel('relative abundance', fontsize = 17)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.legend(fontsize = 14, loc="lower left")
    
    plt.subplot(122)
    plt.scatter(x, y2, label = None)
    plt.scatter(x, y4, label = None)    
    plt.plot(x1, yy2, label = 'Vehicle_2')
    plt.plot(x1, yy4, label = 'Drug_2')
    plt.xlabel('temperture', fontsize = 17)
    plt.ylabel('relative abundance', fontsize = 17)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)  
    plt.legend(fontsize = 14, loc="lower left")


plot_curve(veh_1, cas_1, veh_2, cas_2, 'HSP90AA1')
plot_curve(veh_1, cas_1, veh_2, cas_2, 'HSP90AB1')

