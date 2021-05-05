import numpy as np
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
from tqdm import tqdm
import etienne_tools as et
from astropy.table import Table

import os
from scipy.optimize import curve_fit
import warnings

for ndata in range(1,5):
    alpha = 0.5
    tbl1 = Table.read('data{}_geneva_fp_drift.csv'.format(ndata))
    tbl2 = Table.read('data{}_fp_drift.csv'.format(ndata))

    fig,ax = plt.subplots(nrows = 1, ncols= 1, figsize=[16,8])

    uu = np.unique(tbl1['HIERARCH ESO PRO REC1 CAL9 NAME'])
    for u in uu:
        print(u)
        g = tbl1['HIERARCH ESO PRO REC1 CAL9 NAME'] == u
        plt.vlines(np.nanmax(tbl1['mjdate'][g]), linestyles='--',ymin = -99,ymax = 99,color = 'grey',alpha = 0.5)
        plt.vlines(np.nanmin(tbl1['mjdate'][g]), linestyles='--',ymin = -99,ymax = 99,color = 'grey',alpha = 0.5)

    step = np.nanmedian(np.gradient(np.array(tbl1['mjdate'])))/4.0

    ax.errorbar(tbl1['mjdate']+step,tbl1['rv_A']*1e3,
                 yerr = tbl1['sig_A']*1e3,fmt = '.',color='cyan',label = 'Geneva, A',alpha = alpha)
    ax.errorbar(tbl1['mjdate']+step,tbl1['rv_B']*1e3,
                 yerr = tbl1['sig_B']*1e3,fmt = '.',color='blue',label = 'Geneva, B',alpha = alpha)

    ax.errorbar(tbl2['mjdate'],tbl2['rv_A']*1e3,
                 yerr = tbl2['sig_A']*1e3,fmt = '.',color='red',label = 'UdeM, A',alpha = alpha)
    ax.errorbar(tbl2['mjdate'],tbl2['rv_B']*1e3,
                 yerr = tbl2['sig_B']*1e3,fmt = '.',color='orange',label = 'UdeM, B',alpha = alpha)

    ax.grid(linestyle='--', color='grey', alpha=0.3)

    ax.legend()
    ax.set(xlabel = 'MJDATE',ylabel = 'Velocity [m/s]', title =  ', '.join(np.unique(tbl1['date'])), ylim = [-2,2])
    plt.savefig('dataset_{}.pdf'.format(ndata))
    plt.show()


file = 'data2_fp_drift.csv'
file2 = 'csv_old/data2_fp_drift.csv'
tbl1 = Table.read(file)
tbl2 = Table.read(file2)


