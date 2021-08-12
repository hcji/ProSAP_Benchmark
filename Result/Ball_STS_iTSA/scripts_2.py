# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 09:46:13 2021

@author: jihon
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('Data/Ball_STS_iTSA/Staturosporine_iTSA_data_Ball.csv')
kins = list(data.loc[data.loc[:,'Kinase.Family.Uniprot']=='yes','Accession'].values)

limma_5 = pd.read_csv('Result/Ball_STS_iTSA/limma_52.csv')
limma_4 = pd.read_csv('Result/Ball_STS_iTSA/limma_52_4.csv')
limma_3 = pd.read_csv('Result/Ball_STS_iTSA/limma_52_3.csv')

limma_3_sig = limma_3.loc[limma_3.loc[:, '-logAdjPval'] > 3, 'Accession']
limma_4_sig = limma_4.loc[limma_4.loc[:, '-logAdjPval'] > 3, 'Accession']
limma_5_sig = limma_5.loc[limma_5.loc[:, '-logAdjPval'] > 3, 'Accession']

limma_3_kin = [i for i in limma_3_sig if i in kins]
limma_4_kin = [i for i in limma_4_sig if i in kins]
limma_5_kin = [i for i in limma_5_sig if i in kins]


sig_number = np.array([len(limma_3_sig), len(limma_4_sig), len(limma_5_sig)])
kin_number = np.array([len(limma_3_kin), len(limma_4_kin), len(limma_5_kin)])
non_number = sig_number - kin_number

plot_data = pd.DataFrame({'Method': ['Limma_3_replicates', 'Limma_4_replicates', 'Limma_5_replicates'],
                             'Kinase': kin_number, 'Others': non_number})

plt.figure(figsize=(9, 6.75), dpi = 300)
p1 = plt.bar(np.arange(3), non_number, 0.6, color = '#5B9BD5', bottom = kin_number, edgecolor = 'black')
p2 = plt.bar(np.arange(3), kin_number, 0.6, color = '#C00000', edgecolor = 'black')
plt.tick_params(labelsize = 13)
plt.xticks(np.arange(3), ('Limma_3_replicates', 'Limma_4_replicates', 'Limma_5_replicates'))
plt.ylabel('Numbers', fontsize = 13)
plt.legend((p1[0], p2[0]), ('Others', 'Kinase'), fontsize = 13)
plt.show()


limma_3_ifkin = np.cumsum([i in kins for i in limma_3['Accession']])
limma_4_ifkin = np.cumsum([i in kins for i in limma_4['Accession']])
limma_5_ifkin = np.cumsum([i in kins for i in limma_5['Accession']])

plt.figure(figsize=(6,4.5), dpi=300)
plt.plot(np.arange(1, 101) - limma_3_ifkin[:100], limma_3_ifkin[:100], label='Limma_3_replicates')
plt.plot(np.arange(1, 101) - limma_4_ifkin[:100], limma_4_ifkin[:100], label='Limma_4_replicates')
plt.plot(np.arange(1, 101) - limma_5_ifkin[:100], limma_5_ifkin[:100], label='Limma_5_replicates')
plt.xlabel('Number of false positive')
plt.ylabel('Number of true positive')
plt.legend(loc = 'lower right')

