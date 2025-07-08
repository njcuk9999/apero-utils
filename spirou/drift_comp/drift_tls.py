#!/usr/bin/env python3
#
#  drift.py
#
#  Copyright 2025 operateur <operateur@spip-gis>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
import sys
import os
import optparse
import subprocess
from glob import glob

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.time import Time

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# font = {'size'   : 12}
# plt.rc('font', **font)

# sys.path.append('/mnt/nfs/DRSImages/scripts')
sys.path.append('./')
from read_cfg import read_cfg

locals().update(read_cfg())

outdir = f'{script_out_dir}/drift/'
if not os.path.exists(outdir): os.mkdir(outdir)


########################################################################

def mysavefig(fig, savepath):
    figsize = (6, 3)

    labels_in_fig = True
    for ax in fig.axes:
        formataxis(ax)

    fig.set_size_inches(*figsize)
    fig.set_thick = 100
    fig.tight_layout()

    fig.savefig(savepath + '.pdf', format='pdf')


def formataxis(ax):
    if not ax.get_legend_handles_labels() == ([], []): ax.legend()
    ax.set_xlabel('MJD [days]')
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))


########################################################################

def load_RV(path):
    with fits.open(path) as hdu:
        t = hdu[1].data['MJDMID']
        rv = hdu[1].data['RV'] * 1000
        err = hdu[1].data['DVRMS_SP']
    return t, rv, err


def diff_rv(rv1, err1, rv2, err2):
    rv = rv2 - rv1
    errbar = (err2 ** 2 + err1 ** 2) ** .5
    return rv, errbar


def load_rsb(reddir):
    RSB = {}
    for fiber in ['A', 'B', 'C']:
        files = sorted(glob(f'{reddir}*e2dsff_{fiber}.fits'))
        N = len(files)
        RSB[fiber] = np.zeros(N)
        if fiber == 'C':
            mjd = np.zeros(N)
        for f_idx, file in enumerate(files):
            with fits.open(file) as hdu:
                hdr = hdu[0].header
                RSB[fiber][f_idx] = hdr['EXTSN035f']
                if fiber == 'C':
                    mjd[f_idx] = hdr['MJDMID']
    return mjd, RSB


########################################################################

def load_SS(tstart, tstop, lissage=None):
    """
    tstart : start time in MJD
    tstop : stop time in MJD
    lissage : duration in days for averaging for noise reduction or False for no average
    """
    # tstart = tstart - 3/24

    # ss_tfile = '/mnt/nfs/DRSImages/spip-drs1/FPEnv/FPEnv.txt'
    ss_tfile = './FPEnv.txt'

    t = pd.read_csv(ss_tfile, skiprows=11,
                    names=["time", "fpbody", "fpcase", "fpcase_rms", "tccu", "salle blanche", "controle",
                           "FP pressure (mB)", "Spectro pressure", "Parabola", "Middle bench"])
    t.loc[:, 'time'] = pd.to_datetime(t.loc[:, 'time']).dt.tz_localize("Europe/Paris").dt.tz_convert("UTC")
    mjd = Time(t.loc[:, 'time'].to_list(), scale='utc').mjd  # need to convert to utc
    t.loc[:, 'time'] = mjd

    for col in t: t[col] = pd.to_numeric(t[col], 'coerce')  # convert to float and convert any text to nan

    t = t[(tstart <= t.time) & (t.time < tstop)]
    if not lissage is None:
        t = t.groupby(pd.cut(t['time'], np.arange(tstart, tstop, lissage))).mean()
    return t


########################################################################

def loadRVdict(reddir):
    path = f"{reddir}/cal_drift_FP_FP_<FIBER>.fits"

    # make RV and ERR dicts
    rv_dic = {}
    err_dic = {}
    minmax_dic = {}
    t_dic = {}

    for fiber in ['AB', 'A', 'B', 'C']:
        t, rv, err = load_RV(path.replace('<FIBER>', fiber))
        rv_dic[fiber] = rv
        err_dic[fiber] = err
        t_dic[fiber] = t
        # print("DIC t_dic",fiber,len(t_dic[fiber]))

    # t_AB,rv_AB,err_AB=load_RV(path.replace('<FIBER>','AB')
    # t_A,rv_A,err_A=load_RV(path.replace('<FIBER>','A')
    # t_B,rv_B,err_B=load_RV(path.replace('<FIBER>','B')
    # t_C,rv_C,err_C=load_RV(path.replace('<FIBER>','C')

    # we will go over the fiber C, and check whether the time os there
    # count=0
    # for mjd_C in t_C:
    #    index=np.where
    # for mjd_fiberAB in t_dic['AB']:

    return rv_dic, err_dic, t_dic


def drift(rv_dic, err_dic, t_dic, runname):
    duration = round((t_dic['AB'].max() - t_dic['AB'].min()) * 24, 1)  # h

    errstring = "Fibre" + " |" + "RVstd" + " |" "erbar |" + "MinMax" + "\n"
    errstring = errstring + "------+" * 3 + '------\n'

    fig_rel = plt.figure('relative drift')
    plt.title(f'Relative RV drift over {duration}h')
    ax_rel = plt.subplot(1, 1, 1)
    ax_rel.set_ylabel('RV [m/s]')

    for fiber in ['AB', 'A', 'B']:
        diff, err = diff_rv(rv_dic['C'], err_dic['C'], rv_dic[fiber], err_dic[fiber])
        m_error = np.nanmean(err)
        std_diff = np.nanstd(diff)
        errstring = errstring + f"{fiber + '-C':5s} |{std_diff:.3f} |{m_error:.3f} |{diff.max() - diff.min():.3f}\n"

        ax_rel.plot(t_dic[fiber], diff, label=fiber)

    relpath = f'{outdir}/{runname}_rel_drift'
    mysavefig(fig_rel, relpath)

    fig_abs = plt.figure('absolute drift')
    plt.title(f'Absolute RV drift over {duration}h')
    ax_abs = plt.subplot(1, 1, 1)
    ax_abs.set_ylabel('RV [m/s]')

    for fiber in ['AB', 'A', 'B', 'C']:
        m_error = np.nanmean(err_dic[fiber])
        std_diff = np.nanstd(rv_dic[fiber])
        minmax = rv_dic[fiber].max() - rv_dic[fiber].min()
        errstring = errstring + f"{fiber:5s} |{std_diff:.3f} |{m_error:.3f} |{minmax:.3f}\n"
        ax_abs.plot(t_dic[fiber], rv_dic[fiber], label=fiber)
    abspath = f'{outdir}/{runname}_abs_drift'
    mysavefig(fig_abs, abspath)

    print(errstring)

    outpath = f'{outdir}/{runname}.txt'
    with open(outpath, 'w') as outfile:
        outfile.write(errstring)

    return


########################################################################


def RV_VS_T(rv_dic, err_dic, t_dic, T, runname):
    cols_dic = {"TFP": [('fpcase', 'k')],
                "TSalle": [('salle blanche', 'green'), ('controle', 'purple')],
                "TParab": [('Parabola', 'darkslategrey')],
                "TSpectro": [('Middle bench', 'coral')],
                "TCtrl": [('tccu', 'pink')],
                "P": [('FP pressure (mB)', 'olive')],
                }

    units_dic = {"TFP": 'T [C]',
                 "TSalle": 'T [C]',
                 "TSpectro": 'T [K]',
                 "TParab": 'T [K]',
                 "TCtrl": 'T [C]',
                 "P": 'P [mb]',
                 }

    for figname in ["TFP", "TSalle", "P", "TSpectro", "TParab"]:
        figABS = plt.figure(f'abs vs {figname}')
        axABS = plt.subplot(1, 1, 1)
        plt.title(f'absolute RV vs {figname}')

        figREL = plt.figure(f'rel vs {figname}')
        axREL = plt.subplot(1, 1, 1)
        plt.title(f'relative RV vs {figname}')

        for ax in [axABS, axREL]:
            axT = ax.twinx()
            for col in cols_dic[figname]:  # ('fpcase', 'k'), ('lv1', 'g'),('lv3', 'purple'):
                column = col[0]
                color = col[1]
                axT.plot(T.time, T[column], color=color, linestyle=':', label=column)

            axT.legend()
            axT.set_ylabel(units_dic[figname])

            ax.set_ylabel('RV [m/s]')

        axABS.plot(t_dic['AB'], rv_dic['AB'])
        for i in range(2):  # two empty plots to advance in the color cycle
            axABS.plot([], [])
        axABS.plot(t_dic['C'], rv_dic['C'])
        axREL.plot(t_dic['AB'], rv_dic['AB'] - rv_dic['C'])

        savepath = f'{outdir}/{runname}_drift_vs_{figname}'
        print("Saving :")
        print(savepath + '_abs')
        print(savepath + '_rel')
        mysavefig(figABS, savepath + '_abs')
        mysavefig(figREL, savepath + '_rel')

    return


########################################################################

def main(runname, do_temp=False, do_snr=False, interactive_plot=True):
    reddir = f'{apero_reduced_dir}/{runname}'
    rv_dic, err_dic, t_dic = loadRVdict(reddir)

    ###########################

    drift(rv_dic, err_dic, t_dic, runname)

    plt.close('all')
    if not interactive_plot:
        plt.close('all')

        ###########################
    if do_temp:
        dt = 15 * 60 / 3600 / 24  # s to days # duration on which the temperature is averaged

        tstart = t_dic['AB'].min()  # - 30/3600/24
        tstop = t_dic['AB'].max()  # + 30/3600/24

        T = load_SS(tstart, tstop, lissage=dt)

        RV_VS_T(rv_dic, err_dic, t_dic, T, runname)
    else:
        print('Skipping temperature plot')

    if not interactive_plot:
        plt.close('all')

        ###########################
    if do_snr:
        plt.figure()
        plt.title('SNR')
        plt.xlabel('MJD [days]')
        plt.ylabel('snr')

        mjd, rsb = load_rsb(reddir)
        # print(mjd)
        for lab, val in rsb.items():
            plt.plot(mjd, val, label=lab)
        plt.legend()
        # outdir = '/home/operateur/drs/scripts/out/AT4-10/'
        outdir = './'
        mysavefig(plt.gcf(), f'{outdir}/{runname}_RSB')
    else:
        print('Skipping RSB plot')

    if interactive_plot:
        plt.show()


if __name__ == '__main__':
    parser = optparse.OptionParser()
    # parser.add_option('-a', '--do_abs', default=False, action='store_true')
    parser.add_option('-t', '--do_temp', default=True, action='store_true')
    parser.add_option('-s', '--do_snr', default=False, action='store_true')
    parser.add_option('-i', '--interactive_plot', default=False, action='store_true')
    opt, args = parser.parse_args()

    # do_abs = opt.do_abs
    do_temp = opt.do_temp
    do_snr = opt.do_snr
    interactive_plot = opt.interactive_plot

    runname = args[0]

    main(runname, do_temp, do_snr, interactive_plot)
