#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on {DATE}

@author: cook
"""
import os
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

from apero.core import constants

# =============================================================================
# Define variables
# =============================================================================
__NAME__ = ''
NIGHT = '2020-08-31'

uname = os.environ.get('USERNAME', None)
if uname is None:
    uname = os.environ.get('LOGNAME')

# define where we want to save plots
PLOT_PATH_ALL = dict()
# for rali
if uname == 'spirou':
    PLOT_PATH = '/spirou/cook/paper_plots'
elif uname == 'spirou-client':
    PLOT_PATH = '/home/spirou-client/paper_plots'
elif uname == 'cook':
    PLOT_PATH = '/scratch2/spirou/misc/paper_plots'
else:
    PLOT_PATH = '/data/spirou/drs-data/misc/paper_plots'


# =============================================================================
# Define functions
# =============================================================================
def plot_tcorr_plot(params, case):
    odocode = '2510303'


    if case == 1:
        orders = [33]
        limits_a = [[1665, 1685]]
        limits_b = [[1674, 1676]]
        ylim = [0.6, 1.15]
    else:
        orders = [20]
        limits_a = [[1311, 1316]]
        limits_b = [[1314.5, 1315.5]]
        ylim = [0.4, 1.15]
    # -------------------------------------------------------------------------
    # get filename
    e2dsfile = os.path.join(params['DRS_DATA_OUT'], NIGHT,
                            '{0}e.fits'.format(odocode))
    tcorrfile = os.path.join(params['DRS_DATA_OUT'], NIGHT,
                             '{0}t.fits'.format(odocode))

    # -------------------------------------------------------------------------
    wave = fits.getdata(e2dsfile, extname='WaveAB')
    blaze = fits.getdata(e2dsfile, extname='BlazeAB')
    e2dsffab = fits.getdata(e2dsfile, extname='FluxAB')
    tcorrab = fits.getdata(tcorrfile, extname='FluxAB')
    recon = fits.getdata(tcorrfile, extname='Recon')
    skymodel = fits.getdata(tcorrfile, extname='OHLine')

    # -------------------------------------------------------------------------
    # plot setup
    plt.close()
    plt.figure(figsize=(10, 5))
    shape = (1, 5)
    frame1 = plt.subplot2grid(shape, (0, 0), colspan=3)
    frame2 = plt.subplot2grid(shape, (0, 3), colspan=2)

    frames_right = [frame1]
    frames_left = [frame2]

    # loop around orders
    for it, order_num in enumerate(orders):
        # get this iterations frame
        frame_a = frames_right[it]
        frame_b = frames_left[it]

        frame_part = [frame_a, frame_b]
        limits_part = [limits_a, limits_b]
        # get order values
        ordwave = wave[order_num]
        orde2ds = (e2dsffab / blaze)[order_num]
        ordtcorr = (tcorrab / blaze)[order_num]
        ordrecon = recon[order_num]
        ordskymodel = skymodel[order_num]

        # get first and last point on e2ds
        imask = np.where(np.isfinite(orde2ds))

        # loop around zoom ins
        for part in [0, 1]:

            frame = frame_part[part]
            limits = limits_part[part]

            # get wave limits for first and last poitns of e2ds
            if limits[it] is None:
                wavemin = ordwave[int(np.min(imask))]
                wavemax = ordwave[int(np.max(imask))]
            else:
                wavemin = limits[it][0]
                wavemax = limits[it][1]

            # plot e2dsff
            frame.plot(ordwave, orde2ds / np.nanmedian(orde2ds), color='k',
                       lw=0.5, label='Uncorrected', zorder=2)
            # plot telluric
            frame.plot(ordwave, ordtcorr / np.nanmedian(ordtcorr), color='r',
                       lw=0.5, label='Corrected', zorder=3)
            # plot recon
            frame.plot(ordwave, ordrecon / np.nanmedian(ordrecon), color='b',
                       lw=0.5, label='Telluric absorption', zorder=4,
                       alpha=0.75)

            # set labels
            if part == 0:
                frame.set(xlabel='Normalized flux')
            frame.set(ylim=ylim, xlim=[wavemin, wavemax])
            frame.set(xlabel='Wavelength [nm]')

            # frame.set_title('Order {0}'.format(order_num), y=1.0, pad=-15)

            if it == 0 and part == 0:
                frame.legend(loc='upper center', bbox_to_anchor=(0.85, 1.1),
                             ncol=5)

            # turn off stuff for zoom
            if part == 1:
                waverange = wavemax - wavemin
                point1 = wavemin + 0.25 * waverange
                point2 = wavemax - 0.25 * waverange
                frame.set_xticks([point1, point2])
                frame.set_xticklabels([str(point1), str(point2)])
                frame.set_yticklabels([])

    plt.subplots_adjust(left=0.075, right=0.98, bottom=0.1, top=0.9,
                        hspace=0.25, wspace=0.1)
    # -------------------------------------------------------------------------
    # save plot
    outfile = os.path.join(PLOT_PATH, f'tcorr_farbod_{case}.pdf')
    print('Saving to file: ' + outfile)
    plt.savefig(outfile)
    print('Showing graph')
    plt.show(block=True)
    plt.close()


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # load constants
    _params = constants.load()
    # produce plot
    for _case in [1, 2]:
        plot_tcorr_plot(_params, _case)

# =============================================================================
# End of code
# =============================================================================
