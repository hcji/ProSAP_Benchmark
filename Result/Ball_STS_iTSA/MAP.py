import os
import sys
import time
import base64
import numpy as np
import itertools
import scipy.stats as sps
import matplotlib
matplotlib.rcParams['font.family'] = "Arial"
matplotlib.rc('mathtext', default='regular')
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from statsmodels import api as sm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


Fig_ticks_fontsize = 11
circle_size = 6
line_widths = 0.5

def copyright_show():
    current_date =  time.strftime('%Y', time.localtime())
    implicit = [b'Q29weXJpZ2h0IChjKSAyMDE1LQ==', b'IFNIQU8mWkhBTkcgTGFiLg==', b'QWxsIFJpZ2h0cyBSZXNlcnZlZC4=']
    explicit = [bytes.decode(base64.b64decode(v)) for v in implicit]
    lines = [current_date.join(explicit[:2]), explicit[-1]]
    sys.stdout.write('\n'*2)
    sys.stdout.write('\n'.join(lines))
    sys.stdout.write('\n'*2)
    return None

def progress_bar(pbnum):
    sys.stdout.write('<br>_ %0.1f%% Processed ... s<br><br>' % (100*pbnum))
    return None

def get_data_fun(fn):
    '''
    for get the orignal dataset
    return data_frame, data_header
    '''
    def inner_handle(x):
        return [x[0]] + list(map(float, x[1:3]))
    with open(fn) as cf:
        v = cf.readline().strip().split('\t')
        pre_data = [inner_handle(val.split('\t')) for val in map(str.strip, cf.readlines())]
    hd = v[1:3]
    #sys.stdout.write('\n%s' % v)
    return pre_data, hd

def get_TCnorm_factor_fun(data_sets):
    nfs = []
    now_data = [val[1:] for val in data_sets]
    samp_wise = list(zip(*now_data))
    mean_lib = np.mean(list(map(sum, samp_wise)))
    norm_factor = [mean_lib/sum(val) for val in samp_wise]
    #sys.stdout.write('\n_ #TCNORM %s \n'%len(samp_wise))
    return norm_factor

def get_norm_factor_fun(data_sets, k = 1.5):
    nfs = []
    now_data = [val[1:] for val in data_sets]
    samp_wise = list(zip(*now_data))
    sys.stdout.write('\n_ #%s Samples in Comparison...'%len(samp_wise))
    sys.stdout.write('\n_ #Gene.%s...'%len(samp_wise[0]))

    def plus_one(ele):
        return np.log(ele + 1)
        #return ele
    plot_data = [list(map(plus_one, val)) for val in samp_wise]
    prepare_IQR = [[np.percentile(val, 25), np.percentile(val, 75),
                    np.percentile(val, 75) - np.percentile(val, 25)] for val in plot_data]
    #sys.stdout.write('\n_ IQR')
    #[sys.stdout.write('\n_ %s' % val) for val in prepare_IQR]
    boundarys = [[int(2**(item[0] - k * item[2]) - 1), int(2**(item[1] + (k) * item[2]) - 1)]
                 for item in prepare_IQR]
    #boundarys = [[int(item[0] - k * item[2]), int(item[1] + (k) * item[2])]
    #             for item in prepare_IQR]
    
    valid_data_set = []
    for i, piece in enumerate(now_data):
        flag = list(filter(lambda x: x <0, [ v1 - v2[1] for v1, v2 in zip(piece, boundarys)]))
        if len(flag) == len(piece):
            valid_data_set.append(piece)

    filter_data_sets = list(zip(*valid_data_set))

    sys.stdout.write('<br>_ boundarys')

    lib_size = []
    for idx, sample in enumerate(filter_data_sets):
        lib_size.append(sum(sample))
        sys.stdout.write('<br>_ %s.. ' % boundarys[idx])

        sys.stdout.write('_ %s out of %s ' % (len(sample), len(samp_wise[idx])))
        
    #[sys.stdout.write('%0.3fM'%(ce/1e+6), end = '_| ') for ce in lib_size]
    sys.stdout.write('<br>_ library size... %s' % lib_size)
    
    avg_lib_size = np.mean(lib_size)
    norm_factor = [avg_lib_size/val for val in lib_size]
    sys.stdout.write('<br>_ norm factor.. %s ' % norm_factor)
    return norm_factor

def total_count_normalize_fun(data_sets, norm_factor):
    desire = []
    for index, item in enumerate(data_sets):
        normed_piece = [a*b for a, b in zip(item[1:], norm_factor)]
        desire.append(item[:1] + normed_piece)
    return desire

def get_comp_data(data_sets):
    gene_name, in_data = list(zip(*[[val[0], val[1:]] for val in data_sets]))
    return gene_name, in_data

def get_pairwise_comp_data(data_sets, hd_list, desireidx = [0,-1]):
    cp_title = [ ' vs '.join(val) for idx, val in enumerate(itertools.combinations(hd_list, 2)) if idx in desireidx]
    pairwise_data = [ list(val) for idx, val in enumerate(itertools.combinations(zip(*data_sets), 2)) if idx in desireidx]
    sys.stdout.write('\n_ %s' % cp_title)
    return cp_title, pairwise_data

def get_pairwise_comp_mavalue(data_sets):
    def maget(a,b):
        a = max(a, 1.0)
        b = max(b, 1.0)
        return [np.log2(a/b), 0.5*np.log2(a*b)]
    ma_val = [[ maget(v1, v2) for v1, v2 in zip(*val)] for val in data_sets]
    return ma_val

def get_moving_window_data(ma_data_set, win_size = 500, step_size = 100, center_num = 250, central_site = [0.25, 0.75]):
    # ------------------------win_size filter---------------------------------
    def element_get(ma_paire):
        return '%0.8fININ%0.8f' % (ma_paire[0], ma_paire[1])
    temp_ma_segment_set = set(map(element_get, ma_data_set))
    temp_ma_segment_set_uni = [list(map(float, val.split('ININ'))) for val in temp_ma_segment_set]
       
    sort_ma_set = sorted(filter(lambda x: x[1] != 0, temp_ma_segment_set_uni), key = lambda x: x[1])
    nn = len(sort_ma_set)

    slice_list = []
    for i in range(nn):
        start, end = i * step_size, i * step_size + win_size
        if end < nn:
            slice_list.append([start, end])
        else:
            slice_list.append([nn-win_size, nn])
            break
    sys.stdout.write('\n_ #Ori-%s, #Uni-%s, #sort-%s Length of bins %s, >> %s, %s'%
          (len(ma_data_set), len(temp_ma_segment_set_uni),len(sort_ma_set),len(slice_list),
          slice_list[:2], slice_list[-2:]))
    
    ma_segment_set = [sorted(sort_ma_set[a:b], key = lambda x: x[0]) for a,b in slice_list]
    ma_separate_segment_set = [ list(zip(*val)) for val in ma_segment_set]
    
    format_m_set = [list(val[0]) for val in ma_separate_segment_set]
    a, b = [int(v*win_size) for v in central_site]
    #format_a_set = [np.percentile(val[1][a:b], 10) for val in ma_separate_segment_set]
    format_a_set = [np.median(val[1][a:b]) for val in ma_separate_segment_set]
    return format_m_set, format_a_set

def MS_plotting_position(start_pos, end_pos, total_pos, aa):
    from functools import reduce
    constant = (total_pos - aa + 1) / (total_pos - 2*aa + 1)
    
    def get_inner(x):
        r_digital = (x - aa) / (x - aa + 1)
        return r_digital

    pp_list = []
    for j in range(start_pos, end_pos):
        p_i_list = map(get_inner, range(j, total_pos+1))
        p_i = constant * reduce(lambda m,n: m*n, p_i_list)
        pp_list.append(p_i)
    return pp_list

def quantile_quantile_regression(specific_win_list, percentile_list, theoretical_quantiles, picname = 'null'):
    '''
    return @ variance @ OLS intersept
    '''
    nn = int(len(theoretical_quantiles)/2) + 1
    n1, n2  = nn - 5, nn + 6
    sample_data = sorted(specific_win_list, key = lambda x: x )
    nt = len(sample_data)
    
    #sample_quantiles = [sorted( sample_data[:len(p)] ) for p in theoretical_quantiles]
    sample_quantiles = []
    for p in theoretical_quantiles:
        ccp = len(p)
        a = int(0.5*(nt-ccp))
        b = a + ccp
        sample_quantiles.append(sample_data[a:b])
    
    std_segment, avg_segment = [], []
    std_avg_paired = []
    rsquared_segment = []
    for idix, (tqi, sqj) in enumerate(zip(theoretical_quantiles, sample_quantiles)):
        x = sm.add_constant(tqi)
        y = sqj
        #interceptRLM, slopeRLM = np.nan_to_num(sm.RLM(y, x).fit().params)
        #intercept, slope = sm.OLS(y, x).fit().params
        try:
            results = sm.OLS(y, x).fit()
        except:
            pass
        intercept, slope = list(map(float, np.nan_to_num(results.params)))
        #intercept, slope = np.nan_to_num(sm.RLM(y, x).fit().params)
        rsquared = results.rsquared
        

        std_segment.append(slope)
        avg_segment.append(intercept)
        rsquared_segment.append(rsquared)
        std_avg_paired.append([slope, intercept])
    sorted_std_avg = list(zip(*sorted(std_avg_paired)))
    
    step = min(0.0001, (max(std_segment) - min(std_segment))/26)
    window_list = np.linspace(min(std_segment), max(std_segment)+step, num = 14,endpoint = True)
    moving_list = [[window_list[0], window_list[1]]] + [[item-step, window_list[adx+2]-step] for adx,item in enumerate(window_list[1:-1])]
    #sys.stdout.write(moving_list)
    def overlap_cal(ab, stds):
        a, b = ab
        choose = []
        for v in stds:
            #sys.stdout.write('o',v,end = '')
            if (a-v[0])*(b-v[0]) <= 0:
                choose.append(v)
        #c = input(len(choose))
        return choose
    hist_list= list(map(overlap_cal, moving_list, [std_avg_paired]*len(moving_list) ))
    #sys.stdout.write(hist_list)
    mode_list = []
    for iddx, each in enumerate(hist_list):
        #if iddx % 6 == 0:
        #    sys.stdout.write(len(each),end = '. ')
        if len(each) > len(mode_list):
            mode_list = each
    step_first = list(zip(*mode_list))
    mode = np.median(step_first[0])
    mean = np.mean(step_first[1])

    return (mode**2, mean, rsquared_segment)

def self_mean_var_cal(variance_list_raw):
    variance_list = list(filter(lambda x: abs(x) < 1.3, variance_list_raw))
    n = len(variance_list)
    var = sum(map(np.square, variance_list)) / n
    return var

def exponential_fit_line(xdata, ydata, needx, return_flag='EXP', picshow = False, picname = 'null'):
    from scipy.optimize import curve_fit
    def exp_fun(x, alpha, beta):
        return alpha * np.exp(beta * x)
    #def exp_fun(x, alpha):
    #    return alpha / x

    popt, pcov = curve_fit(exp_fun, xdata, ydata, maxfev=5000)
    #needy = [popt[0] * np.exp(popt[1] * v) for v in  ]    
    needy = popt[0] * np.exp(popt[1]*np.array(needx))
    sys.stdout.write('_ params %s' % popt)
    #needy = [popt[0] / v for v in needx ]

    if picshow:
        plt.plot(xdata, ydata, '>', color = 'indigo', markersize = 15, fillstyle = 'none' )
        plt.plot(needx, needy, '.', color = '#56AEE2', lw = 3, markersize = 3,
                 label = r'$Var$ $=$ $%0.3f\cdot e^{%0.3f \cdot x}$' % (popt[0], popt[1]))
        
        plt.legend(loc = 0, fontsize = 20)
        plt.title(picname)
        plt.ylim(-3,30)
        plt.margins(0.2)
        
        plt.savefig('%s.png'%(picname), dpi=300)
        #plt.show()
        plt.close('all')
 
    if return_flag == 'EXP':
        return needy, popt
    else:
        return list(map(np.log, needy)), (np.log(popt[0]), popt[1])

def significance_get(par_list, data_list):
    def exp_fun(tripar, intensity_raw):
        #c = input(intensity)
        intensity = [val + 1 for val in intensity_raw]
        log_intensity = list(map(np.log2, intensity))
        logfoldchange = log_intensity[0] - log_intensity[1]
        average = np.mean(log_intensity)
        a,b,c = tripar
        std = np.sqrt( a*np.exp(b*average) + c )
        #std = np.sqrt( a*np.exp(b*average) )
        z_score = abs(logfoldchange)/std
        if logfoldchange >= 0:
            p_value = max(1 - sps.norm.cdf(z_score), 1e-50)
        else:
            p_value = max(1 - sps.norm.cdf(z_score) , 1e-50)
            
        return p_value * 2
    p_values = list(map(exp_fun, [par_list]*len(data_list), data_list))
    sys.stdout.write('\n_ testing done...')
    return p_values

def report_result(gene_symbol_list, data_sets, row_pval):
    return_list = []
    def ma_val(signals_raw):
        signals = [val + 1 for val in signals_raw]
        log_signals = list(map(np.log2, signals))
        logfoldchange = log_signals[0] - log_signals[1]
        average = np.mean(log_signals)
        return [logfoldchange, average]
    for idx, item in enumerate(data_sets):
        line = [gene_symbol_list[idx]] + item + ma_val(item) + [row_pval[idx]]
        return_list.append(line)
    return return_list

def report_output(output_name, output_data, output_file = 'output_file'):
    if not os.path.exists(output_file):
        sys.stdout.write( "<br>||%s<br>" % output_file )
        try:
            os.mkdir(output_file)
        except:
            sys.stdout.write( '<br>...%s/%s_MAP.output.xls<br>'%(output_file, output_name) )

    try:
        handle = open('%s/%s_MAP.output.xls'%(output_file, output_name), 'w')
    except:
        sys.stdout.write( 'handle: %s/%s.MAP.output.xls'%(output_file, output_name) )

    a = [handle.write('\t'.join(list(map(str, val))) + '\n') for val in output_data]
    
    handle.close()
    sys.stdout.write('\n_ %s done...\n' % output_name)
    return None

def MA_model_plot(mix_list, model_pars, mean_var_dots, plt_axis, figname = 'null',MAPoutput = 'output_file'):
    if not os.path.exists(MAPoutput):
        os.mkdir(MAPoutput)
    first_mix_list = sorted(mix_list, key = lambda x:x[-1])[::-1]
    merged_list = list(zip(*first_mix_list))
    GS_list, Samp_1, Samp_2, M_value, A_value, P_value = merged_list
    
    mld_y = max(max(M_value), -min(M_value))
    y_ticks = range(int(-mld_y - 1), int(mld_y+3), 2)


    ax = np.linspace(0.0, max(A_value)+1, num = 100)
    ay = model_pars[2] + model_pars[0] * np.exp(ax*model_pars[1])
    vx, vy, vshift = mean_var_dots

    plot_list, model_line = list(zip(*[[i, j] for i, j in zip(ax, ay) if j <= mld_y]))
    plot_list, model_line = np.array(plot_list), np.array(model_line)
                                                                            

    #plot_list = np.linspace(0.5, max(A_value)+1, num = 100)
    #model_line = np.sqrt(model_pars[0]*np.exp(plot_list*model_pars[1]) + model_pars[2]) * abs(sps.norm.ppf(0.05))
    #model_line =  model_pars[2] + model_pars[0] * np.exp(plot_list*model_pars[1])

    colors = [-np.log10( max(val, 1e-5) ) for val in P_value ]
    cfhi = abs(sps.norm.ppf(0.05))
    highlight_idx = [xidx for xidx,iaval in enumerate(A_value) 
                     if abs(M_value[xidx]) / ( np.sqrt(model_pars[0]*np.exp(iaval*model_pars[1]) + model_pars[2]) ) >= cfhi]
    h_xlist = [A_value[hidx] for hidx in highlight_idx]
    h_ylist = [M_value[hidx] for hidx in highlight_idx]

    #=========================================================================================
    plt_axis.plot(A_value, M_value, 'h', markerfacecolor = 'none', markersize = circle_size,
        markeredgewidth = line_widths, markeredgecolor = 'darkgray', label = '%0.5f, %0.3f, %0.3f \n a + b * exp(c * x)' % (model_pars[2], model_pars[0], model_pars[1]))
    plt_axis.plot(h_xlist, h_ylist, '.', markerfacecolor = 'maroon', 
        markeredgewidth = 0.0, markeredgecolor = 'indianred', alpha = 0.8, label = 'Sigs (P-value = 0.05)')
    
    plt_axis.plot(plot_list, model_line, '-', color = 'darkorange', linewidth = line_widths*2, solid_capstyle='round', label = 'Variance trend')
    plt_axis.plot(vx[5:], vy[5:], 'o',markerfacecolor = 'none', markersize = circle_size-3, markeredgewidth = line_widths, 
        markeredgecolor = 'black', label = 'window variance')

    

    plt_axis.plot(vx, vshift, '.-',markerfacecolor = 'none', markersize = circle_size-3, linewidth = line_widths*2, label = 'window shift')

    plt_axis.legend(loc=0, prop={'size':Fig_ticks_fontsize-6})

    plt_axis.plot(plot_list, -model_line, '-', color = 'darkorange', linewidth = line_widths*2,solid_capstyle='round')
    plt_axis.plot([-1, max(A_value)+2], [0,0], '-', color = 'black', linewidth = line_widths,solid_capstyle='round')
    '''
    ax = fig.gca()
    [ax.spines[i].set_linewidth(0.5) for i in ax.spines.keys()]

    ax.yaxis.set_major_locator(MultipleLocator(2.0))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.yaxis.set_minor_locator(MultipleLocator(1.0))
    ax.yaxis.set_tick_params(width = 2, color = 'black', length = 9)
    ax.yaxis.set_tick_params(which ='minor', width = 1.5, length = 5)
    
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_tick_params(width = 2, color = 'black', length = 9, top = 'off')
    ax.xaxis.set_tick_params(which ='minor', width = 1.5, length = 5, top = 'off')
    plt_axis.set_xticks(fontsize = 21, x = -0.3)
    plt_axis.set_yticks( fontsize = 21)
    '''
    x_ticks = range(int(min(A_value)), int(max(A_value)) + 3, 5)
    

    plt_axis.set_xticks(x_ticks)
    plt_axis.set_yticks(y_ticks)
    plt_axis.set_xticklabels(x_ticks, fontdict = {'fontsize': Fig_ticks_fontsize})
    plt_axis.set_yticklabels(y_ticks, fontdict = {'fontsize': Fig_ticks_fontsize})
    
    plt_axis.set_ylabel('Log2 [Fold Changes]', fontsize = Fig_ticks_fontsize, x = -0.05)
    plt_axis.set_xlabel('Mean of log2-intensity', fontsize = Fig_ticks_fontsize, y = -0.05)
    plt_axis.set_title('MA-Plot with the model trend', fontsize = Fig_ticks_fontsize, y = 1.05)

    

    #plt_axis.margins(0.1)
   
    plt_axis.set_ylim(y_ticks[0], y_ticks[-1])
    plt_axis.set_xlim(-1, max(A_value)+2)
    #fig.savefig('%s/%s_MAplot_with_model.pdf'%(MAPoutput,figname))
    #fig.savefig('%s/%s_MAplot_with_model.png'%(MAPoutput,figname))
    sys.stdout.write('<br>_ %s ma model plot done...<br>' % figname)
    #plt.show()
    #plt.close('all')
    return None

def quantile_errer_bar_plot(mlists_of_list, x_quantiles, avg_data, var_data, selected_pect, plt_axis, figname = 'nill', errbar_output =  'output_file', edition = 'Intact'):
    if not os.path.exists(errbar_output):
        os.mkdir(errbar_output)
    def scale_fun(valist, x_bar, x_var):
        return [(val - x_bar)/np.sqrt(x_var) for val in valist]
    
    scaled_data = list(map(scale_fun, mlists_of_list, avg_data, var_data ))
        
    sorted_win_mval = list(map(sorted, scaled_data))
    quantile_wise_mval = list(zip(*sorted_win_mval))
    quantiles_wise_mean = list(map(np.mean, quantile_wise_mval))
    n = len(sorted_win_mval[0])

    # asymmetric_error = [lower_error, upper_error]
    #quantile_errer_bar = list(zip(*[[abs(np.percentile(item, 5) - quantiles_wise_mean[id_2]), abs(quantiles_wise_mean[id_2] - np.percentile(item,95)) ]
    #                                        for id_2, item in enumerate(quantile_wise_mval)]))

    quantile_errer_bar = list(zip(*[[np.std(item, ddof = 1), np.std(item, ddof = 1)] for id_2, item in enumerate(quantile_wise_mval)]))

    annotate_site = sps.norm.ppf(selected_pect)
    #from matplotlib import rcParams
    #rcParams.update({'figure.autolayout': True})
    #import matplotlib.gridspec as gridspec
    start, end = int(n*(1- selected_pect)), int(n*selected_pect)
    
    if edition == 'Intact':
        lx = max(x_quantiles)+0.5
        y_error = quantile_errer_bar

        scatter_size = [20*(a+b) for a,b in zip(*y_error)]
        #f, ax = plt.subplots(121)
        ylim = quantiles_wise_mean[0] - np.std(quantile_wise_mval[0], ddof = 1) - 0.5 , quantiles_wise_mean[-1] + np.std(quantile_wise_mval[-1], ddof = 1)+0.5
        #sys.stdout.write(('quantile_errer_bar[0]',quantile_errer_bar[0]) )

        #plt.subplot(gs[0])
        #plt.subplot(gs[:, :-1])
        #plt.plot([-lx,lx], [0,0],'--', lw = 1, alpha = 0.8, color = 'black', dash_capstyle = 'round')
        #plt.plot([0,0], ylim,'--', lw = 1, alpha = 0.8, color = 'black', dash_capstyle = 'round')
        plt_axis.plot([-lx,lx], [-lx,lx], '-', color = 'orange', lw = line_widths, solid_capstyle = 'round', solid_joinstyle = 'round')

    
        plt_axis.errorbar(x_quantiles, quantiles_wise_mean, yerr=y_error,
                     fmt='.', ecolor = 'royalblue', elinewidth = line_widths*3, capthick = line_widths*3, capsize = line_widths*6, barsabove = False, color = 'dimgrey')
        plt_axis.scatter(x_quantiles, quantiles_wise_mean, s = scatter_size, c = 'dimgrey', marker='o', edgecolors = 'none')

        #plt.plot([annotate_site,annotate_site],[annotate_site-ans_const,annotate_site+ans_const],'--',color = 'dimgray',lw = 4,dash_capstyle = 'round')
        #plt.plot([-annotate_site,-annotate_site],[-annotate_site-ans_const,-annotate_site+ans_const],'--',color = 'dimgray',lw = 4,dash_capstyle = 'round')

        plt_axis.plot([-lx,lx], [-lx,lx], '--', color = 'orange', lw = line_widths*2, dash_capstyle = 'round', dash_joinstyle = 'round')
        plt_axis.plot(x_quantiles[start: end], quantiles_wise_mean[start: end], 'd', color = 'maroon', marker='d', markeredgecolor = 'none')
        
  
        
        plt_axis.set_ylim(ylim[0], ylim[1])

        #leg = plt.legend(loc = [0.40, 0.2], fontsize = 22) # loc = [0.45,0.18]
        #leg.get_frame().set_linewidth(0.0)
        plt_axis.plot([-1,-0.5], [ylim[0]*0.5]*2, 'd', color = 'maroon', markersize = circle_size, markeredgecolor = 'none')
        plt_axis.text(0, ylim[0]*0.53, '%0.1f%% in model' % (200*selected_pect-100.0), fontsize = Fig_ticks_fontsize)


        x_ticks = range(-int(lx), int(lx)+1, 1)
        y_ticks = range(int(ylim[0]), int(ylim[1]), 2)

        plt_axis.set_xticks(x_ticks)
        plt_axis.set_yticks(y_ticks)
        plt_axis.set_xticklabels(x_ticks, fontdict = {'fontsize': Fig_ticks_fontsize})
        plt_axis.set_yticklabels(y_ticks, fontdict = {'fontsize': Fig_ticks_fontsize})

        plt_axis.set_xlim(min(x_quantiles)-0.5, max(x_quantiles)+0.5)
        plt_axis.set_xlabel('Theoretical quantiles',  fontsize = Fig_ticks_fontsize)
        plt_axis.set_ylabel('Unified window quantiles',  fontsize = Fig_ticks_fontsize)
        plt_axis.set_title('\n\n\n%s' % figname, fontsize = Fig_ticks_fontsize, y = 1.05)
    
    #plt.savefig('%s/%s %s.png'%(errbar_output, edition,figname),dpi = 100)
    #progress_bar(1.0)
    #plt.savefig('%s/%s %s.pdf'%(errbar_output, edition,figname))
    #plt.show()
    #plt.close('all')
    sys.stdout.write('\n_ %s ma errbar plot done...\n\n' % figname)
    return None

def MA_pvalue_plot(mix_list, plt_axis, figname = 'null', MAPoutput = 'output_file'):
    if not os.path.exists(MAPoutput):
        os.mkdir(MAPoutput)
    first_mix_list = sorted(mix_list, key = lambda x:x[-1])[::-1]
    merged_list = list(zip(*first_mix_list))
    GS_list, Samp_1, Samp_2, M_value, A_value, P_value = merged_list
    colors = []
    ad_nm = 1e-7
    for val in P_value:
        colors.append(-np.log10(max( val, ad_nm)) )
    
    #fig = plt.figure(figsize = (16,9))   
    psc = plt_axis.scatter(A_value, M_value, s = circle_size*3, c=colors, cmap='OrRd', marker = 'o', edgecolors = 'w',
                facecolors = None, linewidths = line_widths)#, edgecolors = 'face')
    plt_axis.plot([int(min(A_value)) - 1, int(max(A_value))+2], [0, 0], '-', linewidth = 1.0, color = 'black')
    #psc.set_fillstyle = 'none'
    #psc.set_marker = 'h'
    

    '''
    

    ax = fig.gca()
    [ax.spines[i].set_linewidth(2) for i in ax.spines.keys()]

    ax.yaxis.set_major_locator(MultipleLocator(1.0))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax.yaxis.set_tick_params(width = 2, color = 'black', length = 9)
    ax.yaxis.set_tick_params(which ='minor', width = 1.5, length = 5)
        
        #=========================================================================================
    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.xaxis.set_tick_params(width = 2, color = 'black', length = 9, top = 'off')
    ax.xaxis.set_tick_params(which ='minor', width = 1.5, length = 5, top = 'off')
    '''

    x_ticks = range(0, int(max(A_value)) + 3, 5)
    y_ticks = range(int(min(M_value)) - 2, int(max(M_value))+2, 2)

    plt_axis.set_xticks(x_ticks)
    plt_axis.set_yticks(y_ticks)
    plt_axis.set_xticklabels(x_ticks, fontdict = {'fontsize': Fig_ticks_fontsize})
    plt_axis.set_yticklabels(y_ticks, fontdict = {'fontsize': Fig_ticks_fontsize})
    
    plt_axis.set_ylabel('Log2 [Fold Changes]', fontsize = Fig_ticks_fontsize, x = -0.05)
    plt_axis.set_xlabel('Mean of log2-intensity', fontsize = Fig_ticks_fontsize, y = -0.05)
    plt_axis.set_title('MA plot with (-Log10 P-values)', fontsize = Fig_ticks_fontsize, y = 1.05)

    #fig.savefig('%s/%s_MAplot_with_pvlues.pdf'%(MAPoutput,figname))
    #fig.savefig('%s/%s_MAplot_with_pvlues.png'%(MAPoutput,figname))
    #plt.show()
    plt_axis.margins(0.15)
    plt_axis.set_xlim([int(min(A_value)) - 1, int(max(A_value))+2])

    import matplotlib as mpl
    norm = mpl.colors.Normalize(vmin=0, vmax=-np.log10(ad_nm))
    #cb1 = mpl.colorbar.ColorbarBase(plt_axis, cmap='OrRd', norm=norm, orientation='horizontal')
    '''
    plt_axis.axis('off')
    pcol = plt.colorbar(pad = 0.05, fraction = 0.15, format ='%.1f')
    pcol.set_ticks(range(9),['%0.1f'%i for i in range(9)])
    pcol.ax.tick_params(width = 2, length = 9)

    for t in pcol.ax.get_yticklabels():
        t.set_fontsize(18) 
    '''
    

    sys.stdout.write('<br>_ %s ma pvalue plot done...<br>' % figname)
    #plt.close('all')
    return None

if __name__=='__main__':
    sys.stdout.write("<br><br>")
    a = copyright_show()
    sys.stdout.write("<br><br>")

    # -------------------  parameters ---------------------------
    #sys.stdout.write("<br> ------%s-------------  parameters ---------------------------"% (sys.argv[1:]))
    #sys.stdout.write(str(sys.argv[1:]))
    # arvg [1]file_name [2]norm_flag [3]cwin [4]stp [5]central_percentile
    
    file_dir = sys.argv[1][:-1]
    file_name = sys.argv[2]
    norm_flag = sys.argv[3]
    cwin = int(sys.argv[4])
    stp = int(sys.argv[5])
    central_percentile = float(sys.argv[6])/100.0
    '''
    file_dir = '.'
    file_name = 'G:/Work_Protein/Multiple_Pvalue/Raw_Data/UD_forCondModel_527Batch.xls'
    norm_flag = 'Trimmed'
    cwin = 400
    stp = 50
    central_percentile = 50/100.0
    '''
    

    sys.stdout.write('<br><br>_ ------------------------Parameters...--------------<br>')

    sys.stdout.write( '_ 1. File_name.. %s<br>_ 2. Norm_flag.. %s<br>_ 3. Wind_size.. %s<br>_ 4. Step_size.. %s<br>_ 5. Centralpc.. %s%%'%(os.path.basename(file_name), norm_flag, cwin, stp, central_percentile*100) ) 
    sys.stdout.write('<br>_ ------------------------Parameters...--------------<br>')

    pconstant = 0.375
    half_cp = central_percentile * 0.5
    lower_quartile, upper_quartile = 0.5 - half_cp, 0.5 + half_cp
    split_num = int((central_percentile - 0.2) * 100 + 1)
    
    whole_theoretical_quantile = list(map(sps.norm.ppf, MS_plotting_position(1, cwin+1, cwin, pconstant)))

    # 1..................
    progress_bar(0.1)
    # -------------------  parameters ---------------------------

    data_frame, data_title = get_data_fun(file_name)

    if norm_flag == 'pass':
        norm_factors = [1.0, 1.0]
    elif norm_flag == 'Total_count':
        norm_factors = get_TCnorm_factor_fun(data_frame)
    else:
        norm_factors = get_norm_factor_fun(data_frame)
    #sys.stdout.write('<br>============%s=============='% norm_factors)


    norm_data_set = total_count_normalize_fun(data_frame, norm_factors)

    gene_symbol, data_set = get_comp_data(norm_data_set)

    comp_title, pairwise_comp_set = get_pairwise_comp_data(data_set, data_title, desireidx = [0])

    single_ma_value_set = get_pairwise_comp_mavalue(pairwise_comp_set)[0]
    # ------------------------ here just one comp ------------------------
    single_mval_list, single_aval_list = zip(*single_ma_value_set)

    # 2..................
    progress_bar(0.2)

    # --------------------------------------------------------------------
    window_wise_Mval, window_wise_sAval = get_moving_window_data(single_ma_value_set,
                                                                 win_size = cwin, step_size = stp,
                                                                 center_num = int(cwin*central_percentile),
                                                                 central_site = [lower_quartile, upper_quartile])
    win_num = len(window_wise_sAval)
    win_process = list(map(int, np.linspace(0, win_num, num = 3, endpoint = True)))[1:]
    
    # ---------------------------- doing window -------------------------------
    center_percentile_wise = list(np.linspace(0.20, central_percentile, num = 10, endpoint=True))
    
    theoretical_percentile = [MS_plotting_position( int((0.5-p*0.5)*cwin) + 1 , int((0.5+p*0.5)*cwin) + 1, cwin , pconstant)
                              for p in center_percentile_wise]
    theoretical_quantile = [[sps.norm.ppf(valin) for valin in item] for item in theoretical_percentile]
    
    Tot_var_list, Tot_avg_list, Rsqured_list= [], [], []
    for idx, each_win in enumerate(window_wise_Mval):
        temp_var, temp_avg, temp_rs = quantile_quantile_regression(each_win, center_percentile_wise, theoretical_quantile)                
        Tot_var_list.append(temp_var)
        Tot_avg_list.append(temp_avg)
        Rsqured_list.append(temp_rs)
                
        if idx in win_process:
            pass
            #sys.stdout.write( '<br>_ [%s%s] %0.2f %%' % ('-'*idx, ' '*(win_num - idx), 100*(float(idx)/win_num)) )
    #sys.stdout.write(  '<br>_ [%s] 100.0%%' % ('-'*win_num) )
    #---------------------------- window over ----------------------------------
    # 3..................
    progress_bar(0.50)
    
    my_mean_var = self_mean_var_cal(Tot_avg_list)

    var_line_plus, var_paras = exponential_fit_line(window_wise_sAval, Tot_var_list, single_aval_list,
                                                    return_flag='EXP')
    paras_list = list(var_paras) + [my_mean_var] #[0.0]

    # -------------------------- testing -----------------------------
    p_vals = significance_get(paras_list, data_set)

    report_complex = report_result(gene_symbol, data_set, p_vals)


    # ------------------- final report ------------------------------
    progress_bar(0.6)
    file_basename = os.path.basename(file_name)[:-4]
    output_file = '%s/%s' % (file_dir, file_basename + '_MAP.out')
    #input('see') os.path.basename(file_name)[:-4]


    

    header_list = ['#Gene_Symbol'] + ['%s.norm'%val for val in data_title] + ['Mvalue', 'Avalue', 'Pvalue']
    output_data_list = [header_list] + report_complex
    
    report_output(file_basename, output_data_list, output_file = output_file)
    # 4..................
    progress_bar(0.7)

    #  --plot start  v v vv v v v v v v v v v vv v v v  vv v v v v v v vv v  v v vv v v  v v v
    text_color = 'maroon'
    fig, axs = plt.subplots(nrows=1, ncols=3, figsize = (11, 4))
    fig.subplots_adjust(hspace=0.7, wspace= 0.7, left=None, bottom=None, right=None, top=None)
    
    suptitle = os.path.basename(file_name)[:-4] + ' Figs'
    fig.suptitle('%s\n%s' % (suptitle, '::'.join(sys.argv[3:])), fontsize=15, color = text_color)#y = 0.97,   

    progress_bar(0.8)
    sys.stdout.write('\n_ model pars %s\n\n' % paras_list)
    progress_bar(0.9)
    #print Tot_avg_list
    MA_model_plot(report_complex, paras_list, [window_wise_sAval, Tot_var_list, Tot_avg_list], axs[0], figname = file_basename,MAPoutput = output_file)   
    quantile_errer_bar_plot(window_wise_Mval, whole_theoretical_quantile, Tot_avg_list, Tot_var_list, central_percentile*0.5+0.5,
                            axs[1], figname = '%s Quantile Errorbar' % (comp_title[0]), errbar_output =  output_file)

    MA_pvalue_plot(report_complex, axs[2], figname = file_basename, MAPoutput = output_file)
    plt.tight_layout()
    fig.savefig('%s/%s.pdf' % (output_file, suptitle))
    plt.show()
    plt.cla()
    plt.clf()
    plt.close('all')

    #  --plot over ^ ^ ^ ^ ^ ^ ^^ ^ ^ ^ ^ ^ ^^ ^ ^ ^ ^^ ^ ^ ^ 
    progress_bar(1.0)
