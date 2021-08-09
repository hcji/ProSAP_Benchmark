# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 10:53:21 2021

@author: jihon
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from MAP import *


norm_flag = 'Trimmed'
cwin = 400
stp = 50
central_percentile = 50/100.0

pconstant = 0.375
half_cp = central_percentile * 0.5
lower_quartile, upper_quartile = 0.5 - half_cp, 0.5 + half_cp
split_num = int((central_percentile - 0.2) * 100 + 1)


def multi_Map(df, accession_column, group1_columns, group2_columns, log2_transformed = False):
    
    outputs = []
    for i in range(len(group1_columns)):
        g1 = group1_columns[i]
        g2 = group2_columns[i]
        n = df[accession_column].values
        
        if log2_transformed:
            x = 2 ** df[g1].values
            y = 2 ** df[g2].values          
        else:
            x = df[g1].values
            y = df[g2].values
        
        k = np.where((x + y) > 0)[0]
        
        data_title = [g1, g2]
        data_frame = [[n[i], x[i], y[i]] for i in k]
        
        norm_factors = get_norm_factor_fun(data_frame)      
        norm_data_set = total_count_normalize_fun(data_frame, norm_factors)
        gene_symbol, data_set = get_comp_data(norm_data_set)
        comp_title, pairwise_comp_set = get_pairwise_comp_data(data_set, data_title, desireidx = [0])
        single_ma_value_set = get_pairwise_comp_mavalue(pairwise_comp_set)[0]

        single_mval_list, single_aval_list = zip(*single_ma_value_set)
        window_wise_Mval, window_wise_sAval = get_moving_window_data(single_ma_value_set,
                                                                     win_size = cwin, step_size = stp,
                                                                     center_num = int(cwin*central_percentile),
                                                                     central_site = [lower_quartile, upper_quartile])

        win_num = len(window_wise_sAval)
        win_process = list(map(int, np.linspace(0, win_num, num = 3, endpoint = True)))[1:]

        center_percentile_wise = list(np.linspace(0.20, central_percentile, num = 10, endpoint=True))
        theoretical_percentile = [MS_plotting_position( int((0.5-p*0.5)*cwin) + 1 , int((0.5+p*0.5)*cwin) + 1, cwin , pconstant) for p in center_percentile_wise]
        theoretical_quantile = [[sps.norm.ppf(valin) for valin in item] for item in theoretical_percentile]
    
        Tot_var_list, Tot_avg_list, Rsqured_list= [], [], []
        for idx, each_win in enumerate(window_wise_Mval):
            temp_var, temp_avg, temp_rs = quantile_quantile_regression(each_win, center_percentile_wise, theoretical_quantile)                
            Tot_var_list.append(temp_var)
            Tot_avg_list.append(temp_avg)
            Rsqured_list.append(temp_rs)
                
            if idx in win_process:
                pass

        my_mean_var = self_mean_var_cal(Tot_avg_list)
        var_line_plus, var_paras = exponential_fit_line(window_wise_sAval, Tot_var_list, single_aval_list,
                                                    return_flag='EXP')
        paras_list = list(var_paras) + [my_mean_var] #[0.0]

        p_vals = significance_get(paras_list, data_set)
        report_complex = report_result(gene_symbol, data_set, p_vals)

        header_list = ['#Gene_Symbol'] + ['%s.norm'%val for val in data_title] + ['Mvalue', 'Avalue', 'Pvalue']
        output_data_list = pd.DataFrame(report_complex)
        output_data_list.columns = header_list
        outputs.append(output_data_list)
    
    
    output = outputs[0]
    for output_1 in outputs[1:]:
        output = output.merge(output_1, how='inner', on='#Gene_Symbol')
    pv_col = ['Pvalue' in str(i) for i in output.columns]
    pv_mean = np.mean(output.loc[:, pv_col], axis = 1)
    output['AvgPvalue'] = pv_mean
    return output



if __name__=='__main__':
    
    df = pd.read_excel('5_FU_PISA.xlsx')
    accession_column = 'Accession'
    group1_columns = ['DMSO_1', 'DMSO_2', 'DMSO_3', 'DMSO_4', 'DMSO_5']
    group2_columns = ['5FU_1', '5FU_2', '5FU_3', '5FU_4', '5FU_5']
    output = multi_Map(df, accession_column, group1_columns, group2_columns)
    
    
    df = pd.read_csv('Staturosporine_iTSA_data_Ball.csv')
    accession_column = 'Accession'
    group1_columns = ['V_log2.i._TMT_1_iTSA52', 'V_log2.i._TMT_3_iTSA52',
                      'V_log2.i._TMT_5_iTSA52', 'V_log2.i._TMT_7_iTSA52',
                      'V_log2.i._TMT_9_iTSA52']
    group2_columns = ['D_log2.i._TMT_2_iTSA52', 'D_log2.i._TMT_4_iTSA52', 
                      'D_log2.i._TMT_6_iTSA52', 'D_log2.i._TMT_8_iTSA52', 
                      'D_log2.i._TMT_10_iTSA52']
    output = multi_Map(df, accession_column, group1_columns, group2_columns, log2_transformed=True)
    
    map_sig = output.loc[output['AvgPvalue'] < 0.05, '#Gene_Symbol']
    
    kins = list(df.loc[df.loc[:,'Kinase.Family.Uniprot']=='yes','Accession'].values)
    map_kin = [i for i in map_sig if i in kins]
    output.to_csv('D:/project/ProSAP_Benchmark/Result/Ball_STS_iTSA/MAP_52.csv')
    
    