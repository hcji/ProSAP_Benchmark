# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 15:12:09 2021

@author: jihon
"""

import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

data = pd.read_csv('Data/Ball_STS_iTSA/Staturosporine_iTSA_data_Ball.csv')
kins = list(data.loc[data.loc[:,'Kinase.Family.Uniprot']=='yes','Accession'].values)

t_test = pd.read_csv('Result/Ball_STS_iTSA/t_test_52.csv')
limma = pd.read_csv('Result/Ball_STS_iTSA/limma_52.csv')
edgeR = pd.read_csv('Result/Ball_STS_iTSA/edgeR_52.csv')
DESeq = pd.read_csv('Result/Ball_STS_iTSA/DESeq2_52.csv')
Map = pd.read_csv('Result/Ball_STS_iTSA/MAP_52.csv')
Map = Map.sort_values('AvgPvalue')
Map = Map.reset_index(drop=True)


t_test_sig = t_test.loc[t_test.loc[:, '-logAdjPval'] > 3, 'Accession']
limma_sig = limma.loc[limma.loc[:, '-logAdjPval'] > 3, 'Accession']
edgeR_sig = edgeR.loc[edgeR.loc[:, '-logAdjPval'] > 3, 'Accession']
DESeq_sig = DESeq.loc[DESeq.loc[:, '-logAdjPval'] > 3, 'Accession']

t_test_kin = [i for i in t_test_sig if i in kins]
limma_kin = [i for i in limma_sig if i in kins]
edgeR_kin = [i for i in edgeR_sig if i in kins]
DESeq_kin = [i for i in DESeq_sig if i in kins]

t_test_ifkin = np.cumsum([i in kins for i in t_test['Accession']])
limma_ifkin = np.cumsum([i in kins for i in limma['Accession']])
edgeR_ifkin = np.cumsum([i in kins for i in edgeR['Accession']])
DESeq_ifkin = np.cumsum([i in kins for i in DESeq['Accession']])

map_sig = Map.loc[Map.loc[:, 'AvgPvalue'] < 0.001, '#Gene_Symbol']
map_kin = [i for i in map_sig if i in kins]
map_ifkin = np.cumsum([i in kins for i in Map['#Gene_Symbol']])


plt.figure(figsize=(6,4.5), dpi=300)
plt.plot(np.arange(1, 101) - t_test_ifkin[:100], t_test_ifkin[:100], label='T-Test')
plt.plot(np.arange(1, 101) - limma_ifkin[:100], limma_ifkin[:100], label='Limma')
plt.plot(np.arange(1, 101)- edgeR_ifkin[:100], edgeR_ifkin[:100], label='edgeR')
plt.plot(np.arange(1, 101)- DESeq_ifkin[:100], DESeq_ifkin[:100], label='DESeq2')
plt.plot(np.arange(1, 101)- map_ifkin[:100], map_ifkin[:100], label='MAP')
plt.xlabel('Number of false positive')
plt.ylabel('Number of true positive')
plt.legend(loc = 'lower right')
plt.savefig('MAP_1.png')




