# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 10:12:02 2021

@author: hcji
"""

import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sys.path.append(r'D:\project\pyvenn')
import venn
from tqdm import tqdm

data = pd.read_csv('Data/Ball_STS_iTSA/Staturosporine_iTSA_data_Ball.csv')
kins = list(data.loc[data.loc[:,'Kinase.Family.Uniprot']=='yes','Accession'].values)

t_test = pd.read_csv('Result/Ball_STS_iTSA/t_test_52.csv')
limma = pd.read_csv('Result/Ball_STS_iTSA/limma_52.csv')
edgeR = pd.read_csv('Result/Ball_STS_iTSA/edgeR_52.csv')
DESeq = pd.read_csv('Result/Ball_STS_iTSA/DESeq2_52.csv')

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

plt.figure(figsize=(6,4.5), dpi=300)
plt.plot(np.arange(1, 101) - t_test_ifkin[:100], t_test_ifkin[:100], label='T-Test')
plt.plot(np.arange(1, 101) - limma_ifkin[:100], limma_ifkin[:100], label='Limma')
plt.plot(np.arange(1, 101)- edgeR_ifkin[:100], edgeR_ifkin[:100], label='edgeR')
plt.plot(np.arange(1, 101)- DESeq_ifkin[:100], DESeq_ifkin[:100], label='DESeq2')
plt.xlabel('Number of false positive')
plt.ylabel('Number of true positive')
plt.legend(loc = 'lower right')


labels = venn.get_labels([t_test_kin, limma_kin, edgeR_kin, DESeq_kin], fill=['number'])
fig, ax = venn.venn4(labels, names = ['t_test', 'limma', 'edgeR', 'DESeq2'], figsize = (8,6), fontsize=14)
fig.show()


sig_number = np.array([len(t_test_sig), len(limma_sig), len(edgeR_sig), len(DESeq_sig)])
kin_number = np.array([len(t_test_kin), len(limma_kin), len(edgeR_kin), len(DESeq_kin)])
non_number = sig_number - kin_number
    
plot_data = pd.DataFrame({'Method': ['t_test', 'limma', 'edgeR', 'DESeq2'],
                             'Kinase': kin_number, 'Others': non_number})
    
plt.figure(figsize=(6, 4.5), dpi = 300)
p1 = plt.bar(np.arange(4), non_number, 0.6, color = '#5B9BD5', bottom = kin_number, edgecolor = 'black')
p2 = plt.bar(np.arange(4), kin_number, 0.6, color = '#C00000', edgecolor = 'black')
plt.tick_params(labelsize = 13)
plt.xticks(np.arange(4), ('t_test', 'limma', 'edgeR', 'DESeq2'))
plt.ylabel('Numbers', fontsize = 13)
plt.legend((p1[0], p2[0]), ('Others', 'Kinase'), fontsize = 13)
plt.show()


prots = t_test['Accession']
ranks = []
for p in tqdm(prots):
    r1 = 1 + list(t_test['Accession']).index(p)
    r2 = 1 + list(limma['Accession']).index(p)
    r3 = 1 + list(edgeR['Accession']).index(p)
    r4 = 1 + list(DESeq['Accession']).index(p)
    r = (r1 + r2 + r3 + r4) / 4
    ranks.append(r)
consesus = pd.DataFrame(zip(prots, ranks))
consesus.columns = ['Accession', 'Ranks']
consesus = consesus.sort_values('Ranks')
consesus_ifkin = np.cumsum([i in kins for i in consesus['Accession']])

plt.figure(figsize=(6,4.5), dpi=300)
plt.plot(np.arange(1, 101) - t_test_ifkin[:100], t_test_ifkin[:100], label='T-Test')
plt.plot(np.arange(1, 101) - limma_ifkin[:100], limma_ifkin[:100], label='Limma')
plt.plot(np.arange(1, 101)- edgeR_ifkin[:100], edgeR_ifkin[:100], label='edgeR')
plt.plot(np.arange(1, 101)- DESeq_ifkin[:100], DESeq_ifkin[:100], label='DESeq2')
plt.plot(np.arange(1, 101)- consesus_ifkin[:100], DESeq_ifkin[:100], label='Consesus')
plt.xlabel('Number of false positive')
plt.ylabel('Number of true positive')
plt.legend(loc = 'lower right')
