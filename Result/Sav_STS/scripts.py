# -*- coding: utf-8 -*-
"""
Created on Fri May 21 15:33:09 2021

@author: hcji
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.append(r'D:\project\pyvenn')
import venn

data = pd.read_csv('D:/project/CETSA_Benchmark/Raw/Staturosporine_TPP_data_Savitski.csv')
kins = data.loc[data.loc[:,'is kinase?']=='YES','Accession'].values

tpp = pd.read_csv('D:/project/CETSA_Benchmark/Result/Sav_STS/TPP.csv')
tpp_nofilter = pd.read_csv('D:/project/CETSA_Benchmark/Result/Sav_STS/TPP_no_filter.csv')
nparc = pd.read_csv('D:/project/CETSA_Benchmark/Result/Sav_STS/NPARC.csv')
chebyshev = pd.read_csv('D:/project/CETSA_Benchmark/Result/Sav_STS/Chebyshev.csv')
euclidean = pd.read_csv('D:/project/CETSA_Benchmark/Result/Sav_STS/Euclidean.csv')
cosine = pd.read_csv('D:/project/CETSA_Benchmark/Result/Sav_STS/Cosine.csv')
cityblock = pd.read_csv('D:/project/CETSA_Benchmark/Result/Sav_STS/Cityblock.csv')

tpp_sig = tpp.loc[tpp.loc[:,'Score'] > 0,'Accession'].values
tpp_nofilter_sig = tpp_nofilter.loc[tpp_nofilter.loc[:,'Score'] >= 2.417,'Accession'].values
nparc_sig = nparc.loc[chebyshev.loc[:,'Score'] > 3,'Accession'].values
chebyshev_sig = chebyshev.loc[chebyshev.loc[:,'Score'] > 3,'Accession'].values
euclidean_sig = euclidean.loc[euclidean.loc[:,'Score'] > 3,'Accession'].values
cosine_sig = cosine.loc[cosine.loc[:,'Score'] > 3,'Accession'].values
cityblock_sig = cityblock.loc[cityblock.loc[:,'Score'] > 3,'Accession'].values

tpp_kin = [i for i in tpp_sig if i in kins]
nparc_kin = [i for i in nparc_sig if i in kins]
chebyshev_kin = [i for i in chebyshev_sig if i in kins]
euclidean_kin = [i for i in euclidean_sig if i in kins]
cosine_kin = [i for i in cosine_sig if i in kins]
cityblock_kin = [i for i in cityblock_sig if i in kins]

tpp_ifkin = np.cumsum([i in kins for i in tpp['Accession']])
nparc_ifkin = np.cumsum([i in kins for i in nparc['Accession']])
chebyshev_ifkin = np.cumsum([i in kins for i in chebyshev['Accession']])
euclidean_ifkin = np.cumsum([i in kins for i in euclidean['Accession']])
cosine_ifkin = np.cumsum([i in kins for i in cosine['Accession']])
cityblock_ifkin = np.cumsum([i in kins for i in cityblock['Accession']])

plt.figure(figsize=(6,4.5), dpi=300)
plt.plot(np.arange(1, 49) - tpp_ifkin[:48], tpp_ifkin[:48], label='TPP')
plt.plot(np.arange(1, 101) - nparc_ifkin[:100], nparc_ifkin[:100], label='Fitness')
plt.plot(np.arange(1, 101)- nparc_ifkin[:100], chebyshev_ifkin[:100], label='Chebyshev')
plt.plot(np.arange(1, 101)- nparc_ifkin[:100], euclidean_ifkin[:100], label='Euclidean')
plt.plot(np.arange(1, 101)- nparc_ifkin[:100], cosine_ifkin[:100], label='Cosine')
plt.plot(np.arange(1, 101)- nparc_ifkin[:100], cityblock_ifkin[:100], label='Cityblock')
plt.xlabel('Number of false positive')
plt.ylabel('Number of true positive')
plt.legend(loc = 'lower right')

labels = venn.get_labels([tpp_sig, tpp_nofilter_sig, kins], fill=['number'])
plt.figure(dpi = 300)
fig, ax = venn.venn3(labels, names = ['TPP', 'TPP without filter', 'Kinase'], figsize = (8,6), fontsize=14)
fig.show()

labels = venn.get_labels([cosine_kin, cityblock_kin, chebyshev_kin, euclidean_kin], fill=['number'])
plt.figure(dpi = 300)
fig, ax = venn.venn4(labels, names = ['Cosine', 'Cityblock', 'Chebyshev', 'Euclidean'], figsize = (8,6), fontsize=14)
fig.show()


labels = venn.get_labels([tpp_kin, nparc_kin, chebyshev_kin, euclidean_kin], fill=['number'])
plt.figure(dpi = 300)
fig, ax = venn.venn4(labels, names = ['TPP', 'Fitness', 'Chebyshev', 'Euclidean'], figsize = (8,6), fontsize=14)
fig.show()


sig_number = np.array([len(tpp_sig), len(nparc_sig), len(chebyshev_sig), len(euclidean_sig), len(cosine_sig), len(cityblock_sig)])
kin_number = np.array([len(tpp_kin), len(nparc_kin), len(chebyshev_kin), len(euclidean_kin), len(cosine_kin), len(cityblock_kin)])
non_number = sig_number - kin_number
    
plt.figure(figsize=(8,6), dpi = 300)
p1 = plt.bar(np.arange(6), non_number, 0.6, color = '#5B9BD5', bottom = kin_number, edgecolor = 'black')
p2 = plt.bar(np.arange(6), kin_number, 0.6, color = '#C00000', edgecolor = 'black')
plt.tick_params(labelsize = 13)
plt.xticks(np.arange(6), ('TPP', 'Fitness', 'Chebyshev', 'Euclidean', 'Cosine', 'Cityblock'))
plt.yticks(np.arange(0, 80, 10))
plt.ylabel('Numbers', fontsize = 13)
plt.legend((p1[0], p2[0]), ('Others', 'Kinase'), fontsize = 13)
plt.show()

