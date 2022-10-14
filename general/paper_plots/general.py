#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plots for the apero drs paper

Created on 2021-08-01

@author: cook
"""
from astropy.io import fits
from astropy.table import Table
from astropy import constants as cc
from astropy.visualization import imshow_norm, ZScaleInterval
from astropy.visualization import LinearStretch, PowerStretch, LogStretch
import glob
import itertools
import numpy as np
import os
# noinspection PyPep8Naming
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from scipy.signal import convolve2d
from typing import Union
import warnings

import matplotlib; matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Polygon
import matplotlib.patheffects as path_effects

from apero import lang
from apero.core import constants
from apero.core.core import drs_log
from apero.core.core import drs_database
from apero.science import preprocessing as prep
from apero.io import drs_image
from apero.core.instruments.spirou import file_definitions
from apero.science.calib import badpix
from apero.core import math as mp

# =============================================================================
# Define variables
# =============================================================================
__NAME__ = 'paper_plots.general.py'
# Get Logging function
WLOG = drs_log.wlog
# Get the text types
textentry = lang.textentry
# Raw prefix
RAW_PREFIX = file_definitions.raw_prefix
# get the object database
# TODO: remove once on v0.7.243+
if hasattr(drs_database, 'AstrometricDatabase'):
    AstrometricDatabase = drs_database.AstrometricDatabase
else:
    # noinspection PyUnresolvedReferences
    AstrometricDatabase = drs_database.ObjectDatabase
# define the night of data we want to use
NIGHT = '2020-08-31'

uname = os.environ.get('USERNAME', None)
if uname is None:
    uname = os.environ.get('LOGNAME')

# define where we want to save plots
if uname == 'spirou':
    PLOT_PATH = '/spirou/cook/paper_plots'
elif uname == 'cook':
    PLOT_PATH = '/scratch2/spirou/misc/paper_plots'
else:
    PLOT_PATH = '/data/spirou/drs-data/misc/paper_plots'
# define Y, J, H, K
BANDS = dict()
# BANDS['Y'] = [944.126, 1108.771]            # UKIRT Y
# BANDS['J'] = [1080.647, 1406.797]           # 2MASS J
# BANDS['H'] = [1478.738, 1823.102]           # 2MSSS H
# BANDS['K'] = [1954.369, 2344.240]           # 2MASS Ks
BANDS['$Y$'] = [9386.00 / 10, 11134.00 / 10]
BANDS['$J$'] = [11535.86 / 10, 13544.22 / 10]  # MKO
BANDS['$H$'] = [14628.97 / 10, 18085.44 / 10]  # MKO
BANDS['$K_{s}$'] = [19577.92 / 10, 23431.05 / 10]  # MKO

# define plots and append those we want
# noinspection PyListCreation
PLOTS = []

# Apero paper
# PLOTS.append('SIZE_GRID')
# PLOTS.append('FIBER_LAYOUT')
# PLOTS.append('RAW_FEATURES')
# PLOTS.append('PP_FEATURES')
# PLOTS.append('PP_CARTOON')
# PLOTS.append('BADMAP')
# PLOTS.append('BACKMAP')
# PLOTS.append('FLATBLAZE')
# PLOTS.append('E2DS')
# PLOTS.append('S1D')
# PLOTS.append('TCORR')
# PLOTS.append('TELLU_COV')
# PLOTS.append('THERM')
# PLOTS.append('LEAK')

# Telluric paper
PLOTS.append('AIRMASS')

# =============================================================================
# PLOT functions
# =============================================================================
# SIZE_GRID
def plot_size_grid(params):
    odocode = '2510376o'
    hashcode = '2510376o'
    # get file paths
    raw_file = os.path.join(params['DRS_DATA_RAW'], NIGHT,
                            '{0}.fits'.format(odocode))
    pp_file = os.path.join(params['DRS_DATA_WORKING'], NIGHT,
                           '{0}_pp.fits'.format(hashcode))
    e2ds_file = os.path.join(params['DRS_DATA_REDUC'], NIGHT,
                             '{0}_pp_e2dsff_AB.fits'.format(hashcode))
    s1d_file = os.path.join(params['DRS_DATA_REDUC'], NIGHT,
                            '{0}_pp_s1d_v_AB.fits'.format(hashcode))
    ts1d_file = os.path.join(params['DRS_DATA_REDUC'], NIGHT,
                             '{0}_pp_s1d_w_tcorr_AB.fits'.format(hashcode))
    # get
    print('Loading raw image')
    raw_image = fits.getdata(raw_file)
    print('Loading pp image')
    pp_image = fits.getdata(pp_file)
    print('Loading E2DS image')
    e2ds_image = fits.getdata(e2ds_file)
    e2ds_hdr = fits.getheader(e2ds_file)
    print('Loading S1D table')
    s1d_table = Table.read(s1d_file)
    print('Loading tcorr S1D table')
    ts1d_table = Table.read(ts1d_file)
    print('Getting wave sol')
    wavemap = fits2wave(e2ds_image, e2ds_hdr)

    # rotation to match HARPS orientation (expected by DRS)
    # image1 = drs_image.rotate_image(raw_image, params['RAW_TO_PP_ROTATION'])
    # flip image
    image2 = drs_image.flip_image(params, pp_image)
    # get resize size
    sargs = dict(xlow=params['IMAGE_X_LOW'], xhigh=params['IMAGE_X_HIGH'],
                 ylow=params['IMAGE_Y_LOW'], yhigh=params['IMAGE_Y_HIGH'])
    # resize flat
    image3 = drs_image.resize(params, image2, **sargs)

    print('Plotting size_grid plot')
    _ = plt.figure(figsize=(12, 14))
    size = (5, 4)
    frame1 = plt.subplot2grid(size, (0, 0), colspan=2, rowspan=2)
    frame2 = plt.subplot2grid(size, (0, 2), colspan=2, rowspan=2)
    frame3 = plt.subplot2grid(size, (2, 0), colspan=4, rowspan=2)
    frame4 = plt.subplot2grid(size, (4, 0), colspan=4, rowspan=1)

    cmap = matplotlib.cm.get_cmap('inferno').copy()
    cmap.set_bad(color='green')
    # -------------------------------------------------------------------------
    # top left raw image
    _ = _norm_image(raw_image, frame1, cmap)
    # add labels
    frame1.set(xlim=(0, 4096), ylim=(0, 4096))
    frame1.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelleft=False,
                       labelbottom=False)
    frame1.set_title('raw (4096x4096)', loc='left',
                     x=0.05, y=0.95, pad=-14,
                     color='black', backgroundcolor='white')
    # -------------------------------------------------------------------------
    # middle right: flipped + resized image
    _ = _norm_image(image3, frame2, cmap)
    # add labels
    frame2.set(xlim=(-4, 4092), ylim=(0 - 250, 4096 - 250))
    frame2.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelleft=False,
                       labelbottom=False)
    frame2.set_title('pre-processed, flipped, resized (3100x4088)', loc='left',
                     x=0.05, y=0.95, pad=-14,
                     color='black', backgroundcolor='white')
    # -------------------------------------------------------------------------
    # bottom: e2ds
    _ = _norm_image(e2ds_image, frame3, cmap)
    # add labels
    frame3.set(xlim=(0, 4088), ylim=(0, 49))
    frame3.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelleft=False,
                       labelbottom=False)
    frame3.set_title('Extracted (E2DS) 49x4088', loc='left',
                     x=0.025, y=0.95, pad=-14,
                     color='black', backgroundcolor='white')

    # add band regions
    for bandname in BANDS:
        # get band name
        band = BANDS[bandname]
        # wave mask
        wavemask = (wavemap > band[0]) & (wavemap < band[1])
        bandpatch, midpoint = poly_region(wavemask)
        frame3.add_patch(bandpatch)
        # plot text
        txt = frame3.text(midpoint[0], midpoint[1], bandname,
                          color='w', zorder=10, fontsize=16, ha='center')

        txt.set_path_effects([path_effects.withStroke(linewidth=1,
                                                      foreground='k')])

    # -------------------------------------------------------------------------
    # get max flux
    maxflux = np.nanmax(s1d_table['flux'])
    # plot spectrum
    frame4.plot(s1d_table['wavelength'], s1d_table['flux'], lw=0.5,
                color='k', zorder=5)
    frame4.plot(ts1d_table['wavelength'], ts1d_table['flux'], lw=0.5,
                color='r', zorder=10)
    # plot bands
    for bandname in BANDS:
        # get band name
        band = BANDS[bandname]
        # plot band
        frame4.fill_betweenx([-0.2 * maxflux, 1.05 * maxflux], band[0], band[1],
                             color='0.5', alpha=0.5, zorder=0)
        txt = frame4.text(np.mean(band), 0.9 * maxflux, bandname,
                          color='w', zorder=10, fontsize=16, ha='center')
        txt.set_path_effects([path_effects.withStroke(linewidth=1,
                                                      foreground='k')])

    frame4.set_title('Extracted 1D (S1D) 1x285377',
                     loc='left', x=0.5, y=0.95, pad=-14,
                     color='black', backgroundcolor='white')
    # plot cosmetics
    frame4.axes.yaxis.set_ticks([])
    frame4.set(xlim=[950, 2500], ylim=[0, 1.05 * maxflux],
               xlabel='Wavelength [nm]')
    frame4.set_yticklabels([])
    # -------------------------------------------------------------------------
    plt.subplots_adjust(wspace=0.05, hspace=0.05,
                        left=0.01, right=0.99, top=0.975, bottom=0.05)
    # save file
    outfile = os.path.join(PLOT_PATH, 'size_grid.pdf')
    print('Saving to file: ' + outfile)
    plt.savefig(outfile, dpi=300)
    print('Showing graph')
    plt.show()
    plt.close()


# FIBER LAYOUT SCHEMATIC
def plot_fiber_layout(params):
    odocode = '2510301o'
    hashcode = '2510301o'

    zoom0area = [1400, 1700, 1100, 1400]
    zoom1area = [1100, 1400, 1400, 1700]

    # get file paths
    raw_file = os.path.join(params['DRS_DATA_RAW'], NIGHT,
                            '{0}.fits'.format(odocode))
    pp_file = os.path.join(params['DRS_DATA_WORKING'], NIGHT,
                           '{0}_pp.fits'.format(hashcode))

    # get
    print('Loading raw image')
    raw_image = fits.getdata(raw_file)
    print('Loading pp image')
    pp_image = fits.getdata(pp_file)
    # rotation to match HARPS orientation (expected by DRS)
    # image1 = drs_image.rotate_image(raw_image, params['RAW_TO_PP_ROTATION'])
    # flip image
    image2 = drs_image.flip_image(params, pp_image)
    # get resize size
    sargs = dict(xlow=params['IMAGE_X_LOW'], xhigh=params['IMAGE_X_HIGH'],
                 ylow=params['IMAGE_Y_LOW'], yhigh=params['IMAGE_Y_HIGH'])
    # resize flat
    image3 = drs_image.resize(params, image2, **sargs)

    print('Plotting size_grid plot')
    plt.close()

    fig, frames = plt.subplots(ncols=2, nrows=1, figsize=(12, 6))
    frame1 = frames[0]
    frame2 = frames[1]

    cmap = matplotlib.cm.get_cmap('inferno').copy()
    cmap.set_bad(color='green')

    # -------------------------------------------------------------------------
    # zoom raw
    zraw_image = raw_image[zoom0area[2]:zoom0area[3], zoom0area[0]:zoom0area[1]]
    # left raw image
    _, _ = _norm_image(zraw_image, frame1, cmap)
    # add labels
    frame1.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelleft=False,
                       labelbottom=False)

    frame1.set_title('raw', loc='left',
                     x=0.05, y=0.95, pad=-14,
                     color='black', backgroundcolor='white')

    frame1.text(77, 140, 'A', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    frame1.text(94, 115, 'B', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    frame1.text(111, 90, 'C', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)

    frame1.text(77 + 64, 140, 'A', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    frame1.text(94 + 64, 115, 'B', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    frame1.text(111 + 64, 90, 'C', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)

    frame1.text(77 + 126, 140, 'A', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    frame1.text(94 + 126, 115, 'B', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    frame1.text(111 + 126, 90, 'C', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    # -------------------------------------------------------------------------
    # zoom pp
    zimage3 = image3[zoom1area[2]:zoom1area[3], zoom1area[0]:zoom1area[1]]
    # middle: pp
    _, _ = _norm_image(zimage3, frame2, cmap)
    # add labels
    frame2.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelleft=False,
                       labelbottom=False)

    frame2.set_title('pre-processed', loc='left',
                     x=0.05, y=0.95, pad=-14,
                     color='black', backgroundcolor='white')

    frame2.text(100, 234, 'A', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    frame2.text(128, 218, 'B', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    frame2.text(156, 204, 'C', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)

    frame2.text(100, 234 - 53, 'A', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    frame2.text(128, 218 - 53, 'B', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    frame2.text(156, 204 - 53, 'C', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)

    frame2.text(100, 234 - 105, 'A', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    frame2.text(128, 218 - 105, 'B', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)
    frame2.text(156, 204 - 105, 'C', color='black', ha='center', va='center',
                backgroundcolor='white', rotation=7)

    # -------------------------------------------------------------------------
    plt.subplots_adjust(wspace=0.05, hspace=0.05,
                        left=0.01, right=0.99, top=0.975, bottom=0.05)
    # save file
    outfile = os.path.join(PLOT_PATH, 'fiber_layout.pdf')
    print('Saving to file: ' + outfile)
    plt.savefig(outfile, dpi=300)
    print('Showing graph')
    plt.show()
    plt.close()


# RAW_FEATURES
def plot_raw_features(params):
    # set file to use
    odocode = '2510301o'
    # get file paths
    raw_file = os.path.join(params['DRS_DATA_RAW'], NIGHT,
                            '{0}.fits'.format(odocode))
    # get raw image
    print('Loading raw image')
    raw_image = fits.getdata(raw_file)
    # plot feature grid
    plot_feature_grid(raw_image)


def plot_pp_features(params):
    # get pseudo constants
    pconst = constants.pload()
    # set file to use
    odocode = '2510301o'
    # get file paths
    raw_file = os.path.join(params['DRS_DATA_RAW'], NIGHT,
                            '{0}.fits'.format(odocode))
    # get raw image
    print('Loading raw images')
    datalist = []
    for ext in [1, 2, 3, 4]:
        datalist.append(fits.getdata(raw_file, ext=ext))
    # get header
    header = fits.getheader(raw_file)
    # get flux image from the data list
    image = datalist[0]
    # get intercept from the data list
    intercept = datalist[1]
    # get error on slope from the data list
    errslope = datalist[2]
    # get frame time
    frame_time = pconst.FRAME_TIME(params, header)
    # get the pixel exposure time from the data list
    inttime = datalist[3] * frame_time
    # Get hot pixels for corruption check
    hotpixels = prep.get_hot_pixels(params)
    # ----------------------------------------------------------------------
    # Check for pixel shift and/or corrupted files
    # ----------------------------------------------------------------------
    # do this iteratively as if there is a shift need to re-workout QC
    for iteration in range(2):
        # get pass condition
        cout = prep.test_for_corrupt_files(params, image, hotpixels)
        shiftdx, shiftdy = int(cout[2]), int(cout[3])
        # use dx/dy to shift the image back to where the engineering flat
        #    is located
        if shiftdx != 0 and shiftdy != 0:
            # log process
            wmsg = textentry('40-010-00013', args=[shiftdx, shiftdy])
            WLOG(params, '', wmsg)
            # roll on the y axis
            image = np.roll(image, [shiftdy], axis=0)
            intercept = np.roll(intercept, [shiftdy], axis=0)
            errslope = np.roll(errslope, [shiftdy], axis=0)
            inttime = np.roll(inttime, [shiftdy], axis=0)
            # roll on the x axis
            image = np.roll(image, [shiftdx], axis=1)
            intercept = np.roll(intercept, [shiftdx], axis=1)
            errslope = np.roll(errslope, [shiftdx], axis=1)
            inttime = np.roll(inttime, [shiftdx], axis=1)
        elif shiftdx != 0:
            # log process
            wmsg = textentry('40-010-00013', args=[shiftdx, shiftdy])
            WLOG(params, '', wmsg)
            # roll on the x axis
            image = np.roll(image, [shiftdx], axis=1)
            intercept = np.roll(intercept, [shiftdx], axis=1)
            errslope = np.roll(errslope, [shiftdx], axis=1)
            inttime = np.roll(inttime, [shiftdx], axis=1)
        elif shiftdy != 0:
            # log process
            wmsg = textentry('40-010-00013', args=[shiftdx, shiftdy])
            WLOG(params, '', wmsg)
            # roll on the y axis
            image = np.roll(image, [shiftdy], axis=0)
            intercept = np.roll(intercept, [shiftdy], axis=0)
            errslope = np.roll(errslope, [shiftdy], axis=0)
            inttime = np.roll(inttime, [shiftdy], axis=0)
        # ------------------------------------------------------------------
        # correct image
        # ------------------------------------------------------------------
        # correct for the top and bottom reference pixels
        WLOG(params, '', textentry('40-010-00003'))
        image = prep.correct_top_bottom(params, image)
        # correct by a median filter from the dark amplifiers
        WLOG(params, '', textentry('40-010-00004'))
        image = prep.median_filter_dark_amps(params, image)
        # correct for the 1/f noise
        WLOG(params, '', textentry('40-010-00005'))
        image = prep.median_one_over_f_noise(params, image)
        # ---------------------------------------------------------------------
        # Correct for cosmic rays before the possible pixel shift
        # ---------------------------------------------------------------------
        # correct the intercept
        WLOG(params, '', textentry('40-010-00021'))
        intercept = prep.intercept_correct(intercept)
        # correct error slope
        WLOG(params, '', textentry('40-010-00022'))
        errslope1 = prep.errslope_correct(errslope)
        # correct cosmic rays
        WLOG(params, '', textentry('40-010-00018'))
        image, cprops = prep.correct_cosmics(params, image, intercept,
                                             errslope1, inttime)

        # ---------------------------------------------------------------------
        # plot
        # ---------------------------------------------------------------------
        # plot feature grid
        plot_feature_grid(image, 'pp_features.pdf')


def plot_pp_cartoon():
    # size of the image (x and y)
    w = 48

    # get the indices for the image
    y, x = np.indices([w, w])
    # set the reference pixels
    refpix = np.zeros([w, w], dtype=bool)
    refpix[0:4] = True
    refpix[:, 0:4] = True
    refpix[-4:] = True
    refpix[:, -4:] = True
    # add some random noise to the image
    image0 = np.random.normal(0, 1, [w, w])
    # make the random noise lower in the reference pixels
    image0[refpix] -= 3

    # add a fake order to the image
    fit = np.polyfit([0, w / 2.0, w], [3 * w / 4, w / 4, 3 * w / 4], 2)
    expo = 10
    banana1 = np.exp(-0.5 * (x - np.polyval(fit, y)) ** expo / 2 ** expo)
    # add a second order to the image
    bananas = banana1 + np.roll(banana1, 10, axis=1)
    bananas += np.roll(banana1, 20, axis=1)
    # add the orders to the image
    image0[~refpix] += bananas[~refpix] * 5
    # add one over f noise
    npix = 10
    nxpixels = np.arange(npix)
    fspline = ius(nxpixels, np.random.normal(0, 5, [npix]))
    onef = fspline((npix - 2) * np.arange(w) / w + 1) / 2
    # add a second one over f noise signal
    npix = 20
    nxpixels = np.arange(npix)
    fspline = ius(nxpixels, np.random.normal(0, 5, [npix]))
    onef = onef + fspline((npix - 2) * np.arange(w) / w + 1) / 2

    # apply noise to each row of pixels
    for i in range(w):
        image0[:, i] += onef

    # apply amplifier noise to column sections of the detector
    for iamp in range(4):
        start = iamp * (w // 4)
        end = (iamp + 1) * (w // 4)
        image0[:, start:end] += np.random.normal() * 5

    # correct the image for amplifier signal
    per_amp_dc = np.zeros(w * 100)
    image1 = np.array(image0)
    for iamp in range(4):
        start = iamp * (w // 4)
        end = (iamp + 1) * (w // 4)
        med_per_amp = np.nanmedian(image0[0:4, start:end])
        image1[:, start:end] -= med_per_amp
        per_amp_dc[start * 100:end * 100] = med_per_amp

    # correct 1/f
    index = [0, 1, 2, 3, w - 4, w - 3, w - 2, w - 1]
    med = np.nanmedian(image1[:, index], axis=1)
    image2 = np.array(image1)
    for i in range(w):
        image2[:, i] -= med

    # -------------------------------------------------------------------------
    # start plot
    # -------------------------------------------------------------------------
    plt.close()
    fig, frames = plt.subplots(ncols=3, nrows=1, figsize=(12, 4))

    dkwargs = dict(size='15%', pad=0.05)
    # -------------------------------------------------------------------------
    # frame 1
    # -------------------------------------------------------------------------
    frame1 = frames[0]
    # add uncorrected image
    frame1.imshow(image0, origin='lower', aspect='auto')
    # plot dotted line for reference pixels
    frame1.plot(np.array([4, w - 4, w - 4, 4, 4]) - .5,
                np.array([4, 4, w - 4, w - 4, 4]) - .5,
                color='white',
                linestyle='--', alpha=0.6)
    frame1.set_xticks([])
    frame1.set_yticks([])
    frame1.set_xticklabels([])
    frame1.set_yticklabels([])

    divider1 = make_axes_locatable(frame1)
    frame1x = divider1.append_axes('bottom', **dkwargs)
    frame1y = divider1.append_axes('right', **dkwargs)

    frame1x.plot(np.arange(len(per_amp_dc)) / 100, per_amp_dc, color='g',
                 label='Amplifier signal')
    frame1x.set_xticks([])
    frame1x.set_xticklabels([])
    frame1x.set_yticks([])
    frame1x.set_yticklabels([])
    frame1x.set(xlim=[0, len(per_amp_dc) /100],
                ylim=[1.5*np.nanmin(per_amp_dc), 1.5*np.nanmax(per_amp_dc)])

    frame1y.plot(med, np.arange(len(med)), label='Measured 1/f',
                 color='b')
    frame1y.plot(onef, np.arange(len(med)), label='Injected 1/f',
                 color='orange')
    frame1y.set_xticks([])
    frame1y.set_xticklabels([])
    frame1y.set_yticks([])
    frame1y.set_yticklabels([])
    frame1y.set(xlim=[1.5*np.nanmin(onef), 1.5*np.nanmax(onef)],
                ylim=[0, len(med)])
    frame1y.yaxis.tick_right()
    frame1y.yaxis.set_label_position("right")
    # frame1y.legend()
    # -------------------------------------------------------------------------
    # frame 2
    # -------------------------------------------------------------------------
    frame2 = frames[1]

    divider2= make_axes_locatable(frame2)
    frame2x = divider2.append_axes('bottom', **dkwargs)
    frame2y = divider2.append_axes('right', **dkwargs)
    frame2x.axis('off')
    frame2y.axis('off')

    # plot image corrected for amplifiers
    frame2.imshow(image1, origin='lower', aspect='auto')
    # plot dotted line for reference pixels
    frame2.plot(np.array([4, w - 4, w - 4, 4, 4]) - .5,
                np.array([4, 4, w - 4, w - 4, 4]) - .5,
                color='white',
                linestyle='--', alpha=0.6)
    # turn off labels and ticks
    frame2.set_xticks([])
    frame2.set_yticks([])
    frame2.set_xticklabels([])
    frame2.set_yticklabels([])
    # -------------------------------------------------------------------------
    # frame 2
    # -------------------------------------------------------------------------
    frame3 = frames[2]

    divider3= make_axes_locatable(frame3)
    frame3x = divider3.append_axes('bottom', **dkwargs)
    frame3y = divider3.append_axes('right', **dkwargs)
    frame3x.axis('off')
    frame3y.axis('off')

    frame3.imshow(image2, origin='lower', aspect='auto')
    # plot dotted line for reference pixels
    frame3.plot(np.array([4, w - 4, w - 4, 4, 4]) - .5,
                np.array([4, 4, w - 4, w - 4, 4]) - .5,
                color='white',
                linestyle='--', alpha=0.6)
    # turn off labels and ticks
    frame3.set_xticks([])
    frame3.set_yticks([])
    frame3.set_xticklabels([])
    frame3.set_yticklabels([])
    # -------------------------------------------------------------------------
    plt.subplots_adjust(left=0.01, right=0.99, bottom=0.05, top=0.975,
                        hspace=0.05, wspace=0.05)

    # add legend
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    fig.legend(lines, labels, loc=(0.35, 0.05), ncol=3)

    # save file
    outfile = os.path.join(PLOT_PATH, 'pp_cartoon.pdf')
    print('Saving to file: ' + outfile)
    plt.savefig(outfile, dpi=300)
    print('Showing graph')
    plt.show()
    plt.close()


# plot for raw features and pp features
def plot_feature_grid(image, outname='raw_features.pdf'):
    # define cuts / zooms
    cut1area = [100, 700, 1200, 1800]
    cut2area = [2375, 2575, 525, 725]
    zoom0area = [800, 1200, 3396, 3796]
    zoom1area = [3596, 3996, 800, 1200]
    zoom2area = [1240, 1290, 1590, 1640]
    # text height above image [in raw image pixel units]
    textheight = 100
    rcolor = 'r'
    # -------------------------------------------------------------------------
    print('Plotting raw_features plot')
    # set up grid
    plt.close()
    fig = plt.figure(figsize=(12, 14))
    size = (3, 3)
    topright = plt.subplot2grid(size, (0, 0))
    topmid = plt.subplot2grid(size, (0, 1))
    topleft = plt.subplot2grid(size, (0, 2))
    rightmid = plt.subplot2grid(size, (1, 0))
    rightbot = plt.subplot2grid(size, (2, 0))
    panel = plt.subplot2grid(size, (1, 1), colspan=2, rowspan=2)

    # get colour map
    cmap = matplotlib.cm.get_cmap('inferno').copy()
    cmap.set_bad(color='green')

    cmap0 = matplotlib.cm.get_cmap('Greys_r').copy()
    cmap0.set_bad(color='green')
    # -------------------------------------------------------------------------
    # dark region
    cut1 = image[cut1area[2]:cut1area[3], cut1area[0]:cut1area[1]]
    c1im, _ = _norm_image(cut1, rightmid, cmap)
    _add_colorbar(fig, c1im, rightmid, side='bottom')
    rightmid.set_title('D. dark region', loc='left',
                       x=0.05, y=0.95, pad=-14,
                       color='black', backgroundcolor='white')
    rightmid.tick_params(axis='both', which='both', bottom=False, top=False,
                         left=False, right=False, labelleft=False,
                         labelbottom=False)
    # -------------------------------------------------------------------------
    # holes
    cut2 = image[cut2area[2]:cut2area[3], cut2area[0]:cut2area[1]]
    c2im, _ = _norm_image(cut2, rightbot, cmap)
    _add_colorbar(fig, c2im, rightbot, side='bottom')
    rightbot.set_title('E. detector holes', loc='left',
                       x=0.05, y=0.95, pad=-14,
                       color='black', backgroundcolor='white')
    rightbot.tick_params(axis='both', which='both', bottom=False, top=False,
                         left=False, right=False, labelleft=False,
                         labelbottom=False)
    # -------------------------------------------------------------------------
    # zoom 0
    zoom0 = image[zoom0area[2]:zoom0area[3], zoom0area[0]:zoom0area[1]]
    z0im, _ = _norm_image(zoom0, topright, cmap)
    _add_colorbar(fig, z0im, topright, side='bottom')
    topright.set_title('A. Zoom reddest orders', loc='left',
                       x=0.05, y=0.95, pad=-14,
                       color='black', backgroundcolor='white')
    topright.tick_params(axis='both', which='both', bottom=False, top=False,
                         left=False, right=False, labelleft=False,
                         labelbottom=False)
    # -------------------------------------------------------------------------
    # zoom 1
    zoom1 = image[zoom1area[2]:zoom1area[3], zoom1area[0]:zoom1area[1]]
    z1im, norm = _norm_image(zoom1, topmid, cmap)
    _add_colorbar(fig, z1im, topmid, side='bottom')
    topmid.set_title('B. Zoom bluest orders', loc='left',
                     x=0.05, y=0.95, pad=-14,
                     color='black', backgroundcolor='white')
    topmid.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelleft=False,
                       labelbottom=False)
    # -------------------------------------------------------------------------
    # zoom 2
    zoom2 = image[zoom2area[2]:zoom2area[3], zoom2area[0]:zoom2area[1]]
    z2im, norm = _norm_image(zoom2, topleft, cmap)
    _add_colorbar(fig, z2im, topleft, side='bottom')
    topleft.set_title('C. Zoom on slices', loc='left',
                      x=0.05, y=0.95, pad=-14,
                      color='black', backgroundcolor='white')
    topleft.tick_params(axis='both', which='both', bottom=False, top=False,
                        left=False, right=False, labelleft=False,
                        labelbottom=False)
    # -------------------------------------------------------------------------
    # panel
    pim, norm = _norm_image(image, panel, cmap0)
    _add_colorbar(fig, pim, panel, side='bottom')
    panel.tick_params(axis='both', which='both', bottom=False, top=False,
                      left=False, right=False, labelleft=False,
                      labelbottom=False)
    panel.set_title('raw (4096x4096)', loc='left',
                    x=0.5, y=0.95, pad=-14,
                    color='black', backgroundcolor='white')
    # -------------------------------------------------------------------------
    # add rectangle dark cut
    rcut1 = patches.Rectangle((cut1area[0], cut1area[2]),
                              cut1area[1] - cut1area[0],
                              cut1area[3] - cut1area[2],
                              linewidth=1, edgecolor=rcolor, facecolor='none')
    panel.add_patch(rcut1)
    # add text
    panel.text(0.5 * (cut1area[1] + cut1area[0]), cut1area[3] + textheight,
               'D', color=rcolor, backgroundcolor='white')
    # -------------------------------------------------------------------------
    # add rectangle holes
    rcut2 = patches.Rectangle((cut2area[0], cut2area[2]),
                              cut2area[1] - cut2area[0],
                              cut2area[3] - cut2area[2],
                              linewidth=1, edgecolor=rcolor, facecolor='none')
    panel.add_patch(rcut2)
    # add text
    panel.text(0.5 * (cut2area[1] + cut2area[0]), cut2area[3] + textheight,
               'E', color=rcolor, backgroundcolor='white')
    # -------------------------------------------------------------------------
    # add rectangle zoom 1
    rzoom0 = patches.Rectangle((zoom0area[0], zoom0area[2]),
                               zoom0area[1] - zoom0area[0],
                               zoom0area[3] - zoom0area[2],
                               linewidth=1, edgecolor=rcolor, facecolor='none')
    panel.add_patch(rzoom0)
    # add text
    panel.text(0.5 * (zoom0area[1] + zoom0area[0]), zoom0area[3] + textheight,
               'A', color=rcolor, backgroundcolor='white')
    # -------------------------------------------------------------------------
    # add rectangle zoom 1
    rzoom1 = patches.Rectangle((zoom1area[0], zoom1area[2]),
                               zoom1area[1] - zoom1area[0],
                               zoom1area[3] - zoom1area[2],
                               linewidth=1, edgecolor=rcolor, facecolor='none')
    panel.add_patch(rzoom1)
    # add text
    panel.text(0.5 * (zoom1area[1] + zoom1area[0]), zoom1area[3] + textheight,
               'B', color=rcolor, backgroundcolor='white')
    # -------------------------------------------------------------------------
    # add rectangle zoom 2
    rzoom2 = patches.Rectangle((zoom2area[0], zoom2area[2]),
                               zoom2area[1] - zoom2area[0],
                               zoom2area[3] - zoom2area[2],
                               linewidth=1, edgecolor=rcolor, facecolor='none')
    panel.add_patch(rzoom2)
    # add text
    panel.text(0.5 * (zoom2area[1] + zoom2area[0]), zoom2area[3] + textheight,
               'C', color=rcolor, backgroundcolor='white')
    # -------------------------------------------------------------------------
    plt.subplots_adjust(wspace=0.05, hspace=0.075,
                        left=0.01, right=0.99, top=0.975, bottom=0.025)
    # -------------------------------------------------------------------------
    outfile = os.path.join(PLOT_PATH, outname)
    print('Saving to file: ' + outfile)
    plt.savefig(outfile, dpi=300)
    print('Showing graph')
    plt.show()
    plt.close()


# BADMAP
def plot_badpix_plot(params):
    hashcode = 'DB67D5C4F5'
    xlow = params['IMAGE_X_LOW']
    xhigh = params['IMAGE_X_HIGH']
    ylow = params['IMAGE_Y_LOW']
    yhigh = params['IMAGE_Y_HIGH']

    dark_file = glob.glob(os.path.join(params['DRS_DATA_REDUC'], 'other',
                                       '*d_pp_dark_master.fits'))[0]
    bad_file = os.path.join(params['DRS_DATA_REDUC'], NIGHT,
                            '{0}_pp_badpixel.fits'.format(hashcode))

    dark_image = fits.getdata(dark_file)
    bad_image = fits.getdata(bad_file).astype(float)

    # fill bad image
    bad_image_full = np.zeros_like(dark_image)
    bad_image_full[ylow:yhigh, xlow:xhigh] = bad_image
    dark_image = drs_image.flip_image(params, dark_image)

    cmap1 = matplotlib.cm.get_cmap('Greys_r').copy()
    cmap2 = matplotlib.cm.get_cmap('Greys').copy()

    plt.close()
    fig, frames = plt.subplots(figsize=(20, 10), ncols=2, nrows=1)
    frame0 = frames[0]
    frame1 = frames[1]

    _, _ = imshow_norm(dark_image, frame0, origin='lower', aspect='auto',
                       interval=ZScaleInterval(), stretch=LinearStretch(),
                       cmap=cmap1, interpolation='None', rasterized=True)
    _, _ = imshow_norm(bad_image_full, frame1, origin='lower', aspect='auto',
                       cmap=cmap2, interpolation='None', rasterized=True)

    frame0.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelleft=False,
                       labelbottom=False)
    frame1.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelleft=False,
                       labelbottom=False)

    title1 = ('Original shape: ({0}x{1})  {2:,d} total pixels'
              '\nResized shape: ({3}x{4})  {5:,d} total pixels')
    title2 = ('Number of bad pixels (in resized region): = {0:,d}'
              '\nFraction of bad pixels (in resized region): = {1:.2f} %')

    farg1 = list(dark_image.shape) + [np.product(dark_image.shape)]
    farg1 += list(bad_image.shape) + [np.product(bad_image.shape)]
    frame0.set_title(title1.format(*farg1), y=1.0, pad=-100,
                     backgroundcolor='white', fontsize=16)
    farg2 = [int(np.sum(bad_image_full)),
             100 * np.sum(bad_image_full) / np.product(bad_image.shape)]
    frame1.set_title(title2.format(*farg2), y=1.0, pad=-100,
                     backgroundcolor='silver', fontsize=16)

    frame0.hlines(y=ylow, xmin=xlow, xmax=xhigh, color='r', lw=2)
    frame0.hlines(y=yhigh, xmin=xlow, xmax=xhigh, color='r', lw=2)
    frame0.vlines(x=xlow, ymin=ylow, ymax=yhigh, color='r', lw=2)
    frame0.vlines(x=xhigh, ymin=ylow, ymax=yhigh, color='r', lw=2)
    frame1.hlines(y=ylow, xmin=xlow, xmax=xhigh, color='r', lw=2)
    frame1.hlines(y=yhigh, xmin=xlow, xmax=xhigh, color='r', lw=2)
    frame1.vlines(x=xlow, ymin=ylow, ymax=yhigh, color='r', lw=2)
    frame1.vlines(x=xhigh, ymin=ylow, ymax=yhigh, color='r', lw=2)

    frame0.set(xlim=[0, 4096], ylim=[0, 4096])
    frame1.set(xlim=[0, 4096], ylim=[0, 4096])

    plt.subplots_adjust(wspace=0, hspace=0, left=0.01, right=0.99,
                        bottom=0.01, top=0.99)

    outfile = os.path.join(PLOT_PATH, 'badmap.pdf')
    print('Saving to file: ' + outfile)
    plt.savefig(outfile)
    print('Showing graph')
    plt.show()
    plt.close()


# BACKMAP
def plot_backmap_plot(params):
    # get constants
    width = params['BKGR_BOXSIZE']
    # set a pid
    params['PID'] = 'TEST'
    # define hash codes
    flathash = 'DB67D5C4F5'
    darkhash = 'C3EB202753'
    # -------------------------------------------------------------------------
    # get filenames
    flat_file = os.path.join(params['DRS_DATA_REDUC'], NIGHT,
                             '{0}_pp.fits'.format(flathash))
    dark_file = os.path.join(params['DRS_DATA_REDUC'], NIGHT,
                             '{0}_pp.fits'.format(darkhash))
    # get images
    flat_image = fits.getdata(flat_file)
    dark_image = fits.getdata(dark_file)
    # get shape of original image
    nshape = flat_image.shape
    # -------------------------------------------------------------------------
    # run the bad pixel code relevant for the plot
    WLOG(params, 'info', 'Running badpix code required for plot')
    flat_image1, backmask, backest = _do_bad_pix(params, flat_image, dark_image)

    # set 0.0 back to NaN for graph
    flat_image1[flat_image1 == 0.0] = np.nan
    # -------------------------------------------------------------------------
    # get colour maps
    cmap1 = matplotlib.cm.get_cmap('inferno').copy()
    cmap2 = matplotlib.cm.get_cmap('Greys_r').copy()
    cmap3 = matplotlib.cm.get_cmap('Greys').copy()
    # set NaNs to green
    cmap2.set_bad(color='green')
    # -------------------------------------------------------------------------
    # plot setup
    plt.close()
    fig, frames = plt.subplots(figsize=(18, 6), ncols=3, nrows=1)
    frame0 = frames[0]
    frame1 = frames[1]
    frame2 = frames[2]
    # -------------------------------------------------------------------------
    # three imshow norm plots
    _, _ = imshow_norm(flat_image1, frame0, origin='lower', aspect='auto',
                       interval=ZScaleInterval(), stretch=LinearStretch(),
                       cmap=cmap2, interpolation='None', rasterized=True)
    _, _ = imshow_norm(backest, frame1, origin='lower', aspect='auto',
                       cmap=cmap1, interpolation='None', rasterized=True)
    _, _ = imshow_norm(backmask, frame2, origin='lower', aspect='auto',
                       cmap=cmap3, interpolation='None', rasterized=True)
    # -------------------------------------------------------------------------
    # plot vertical lines on frame0
    lines = np.arange(0, nshape[0] + width, width) - 4
    # remove first and last
    lines = lines[1:-1]
    for line in lines:
        frame0.axvline(line, color='r')
    # -------------------------------------------------------------------------
    # remove ticks
    frame0.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelleft=False,
                       labelbottom=False)
    frame1.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelleft=False,
                       labelbottom=False)
    frame2.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelleft=False,
                       labelbottom=False)

    frame2.set(xlim=[0, 512], ylim=[512, 1024])

    # adjust edges of figure
    plt.subplots_adjust(wspace=0, hspace=0, left=0.01, right=0.99,
                        bottom=0.01, top=0.99)
    # -------------------------------------------------------------------------
    # save plot
    outfile = os.path.join(PLOT_PATH, 'backmap.pdf')
    print('Saving to file: ' + outfile)
    plt.savefig(outfile)
    print('Showing graph')
    plt.show()
    plt.close()


# FLAT BLAZE
def plot_flatblaze_plot(params):
    # set a pid
    params['PID'] = 'TEST'
    # define hash codes
    flathash = 'DB67D5C4F5'
    # list orders to plot
    orders = [5, 15, 25, 35, 45]

    colors = ['orange', 'r', 'g', 'b', 'purple']

    # -------------------------------------------------------------------------
    # get filenames
    flat_file = os.path.join(params['DRS_DATA_REDUC'], NIGHT,
                             '{0}_pp_flat_AB.fits'.format(flathash))
    blaze_file = os.path.join(params['DRS_DATA_REDUC'], NIGHT,
                              '{0}_pp_blaze_AB.fits'.format(flathash))
    # get images
    flat_image = fits.getdata(flat_file)
    blaze_image = fits.getdata(blaze_file)
    e2ds_image = flat_image * blaze_image
    # -------------------------------------------------------------------------
    # plot setup
    plt.close()
    _ = plt.figure(figsize=(14, 8))

    frame0 = plt.subplot2grid((3, 1), (0, 0), rowspan=2, colspan=1)
    frame1 = plt.subplot2grid((3, 1), (2, 0), rowspan=1, colspan=1)

    for it, order_num in enumerate(orders):
        order_name = 79 - order_num

        # norm = np.nanmedian(e2ds_image[order_num])

        frame0.plot(e2ds_image[order_num], label=f'Order #{order_name}',
                    color=colors[it], alpha=0.5)
        frame0.plot(blaze_image[order_num], color=colors[it], linestyle='--')
        frame1.plot(flat_image[order_num], color=colors[it], alpha=0.75)

    frame0.legend(loc=0)
    frame0.set(xlabel='Pixel number', ylabel='E2DS Flux [e-]',
               xlim=[0, 4088])
    frame1.set(xlabel='Pixel number', ylabel='Flat = E2DS/Blaze',
               xlim=[0, 4088])
    # adjust edges of figure
    plt.subplots_adjust(wspace=0, hspace=0.225, left=0.075, right=0.99,
                        bottom=0.075, top=0.95)
    # -------------------------------------------------------------------------
    # save plot
    outfile = os.path.join(PLOT_PATH, 'flat_blaze.pdf')
    print('Saving to file: ' + outfile)
    plt.savefig(outfile)
    print('Showing graph')
    plt.show(block=True)
    plt.close()


# E2DS
def plot_e2ds_plot(params):
    odocode = '2510303o'
    # -------------------------------------------------------------------------
    # get filenames
    pp_file = os.path.join(params['DRS_DATA_WORKING'], NIGHT,
                           '{0}_pp.fits'.format(odocode))
    e2ds_file = os.path.join(params['DRS_DATA_REDUC'], NIGHT,
                             '{0}_pp_e2ds_AB.fits'.format(odocode))
    e2dsll_file = os.path.join(params['DRS_DATA_REDUC'], NIGHT,
                               'DEBUG_{0}_pp_e2dsll_AB.fits'.format(odocode))
    # -------------------------------------------------------------------------
    # get images
    pp_image = fits.getdata(pp_file)
    # get resize size
    sargs = dict(xlow=params['IMAGE_X_LOW'],
                 xhigh=params['IMAGE_X_HIGH'],
                 ylow=params['IMAGE_Y_LOW'],
                 yhigh=params['IMAGE_Y_HIGH'])
    # flip image
    pp_image = drs_image.flip_image(params, pp_image)
    # resize pp image
    pp_image = drs_image.resize(params, pp_image, **sargs)

    e2dsll_image = fits.getdata(e2dsll_file)
    e2ds_image = fits.getdata(e2ds_file)

    # get zoom
    ylow0, yhigh0 = 1450, 1650
    xlow0, xhigh0 = 700, 850

    pp_image1 = pp_image[ylow0:yhigh0, xlow0:xhigh0]

    # scale box to e2dsll
    ylow1a, yhigh1a = 825, 925
    e2dsll_image1 = e2dsll_image[ylow1a:yhigh1a, xlow0:xhigh0]

    # scale box to e2dsll
    ylow1b, yhigh1b = 25, 30

    e2ds_image1 = e2ds_image[ylow1b:yhigh1b, xlow0:xhigh0]

    # -------------------------------------------------------------------------
    # get colour maps
    cmap1 = matplotlib.cm.get_cmap('inferno').copy()
    cmap1.set_bad(color='green')
    # -------------------------------------------------------------------------
    # plot setup
    plt.close()
    fig, frames = plt.subplots(figsize=(18, 12), ncols=3, nrows=2)
    frame0 = frames[0][0]
    frame1 = frames[0][1]
    frame2 = frames[0][2]
    frame3 = frames[1][0]
    frame4 = frames[1][1]
    frame5 = frames[1][2]
    # -------------------------------------------------------------------------
    # three imshow norm plots
    _, _ = imshow_norm(pp_image, frame0, origin='lower', aspect='auto',
                       interval=ZScaleInterval(), stretch=LinearStretch(),
                       cmap=cmap1, interpolation='None', rasterized=True)
    _, _ = imshow_norm(e2dsll_image, frame1, origin='lower', aspect='auto',
                       cmap=cmap1, interpolation='None', rasterized=True)
    _, _ = imshow_norm(e2ds_image, frame2, origin='lower', aspect='auto',
                       cmap=cmap1, interpolation='None', rasterized=True)

    _, _ = imshow_norm(pp_image1, frame3, origin='lower', aspect='auto',
                       interval=ZScaleInterval(), stretch=LinearStretch(),
                       cmap=cmap1, interpolation='None', rasterized=True)
    _, _ = imshow_norm(e2dsll_image1, frame4, origin='lower', aspect='auto',
                       cmap=cmap1, interpolation='None', rasterized=True)
    _, _ = imshow_norm(e2ds_image1, frame5, origin='lower', aspect='auto',
                       cmap=cmap1, interpolation='None', rasterized=True)
    # -------------------------------------------------------------------------
    # remove ticks
    for frame in frames.ravel():
        frame.tick_params(axis='both', which='both', bottom=False, top=False,
                          left=False, right=False, labelleft=False,
                          labelbottom=False)
    # adjust edges of figure
    plt.subplots_adjust(wspace=0, hspace=0, left=0.01, right=0.99,
                        bottom=0.01, top=0.99)
    # -------------------------------------------------------------------------
    # save plot
    outfile = os.path.join(PLOT_PATH, 'e2ds_grid.pdf')
    print('Saving to file: ' + outfile)
    plt.savefig(outfile)
    print('Showing graph')
    plt.show()
    plt.close()


# S1D
def plot_s1d_plot(params):
    odocode = '2510303o'
    first_order = 34
    last_order = 36
    # -------------------------------------------------------------------------
    # get filenames
    e2ds_file = os.path.join(params['DRS_DATA_REDUC'], NIGHT,
                             '{0}_pp_e2dsff_AB.fits'.format(odocode))
    # -------------------------------------------------------------------------
    # get images
    sp1 = fits.getdata(e2ds_file)
    hdr = fits.getheader(e2ds_file)
    # -------------------------------------------------------------------------
    # get blaze file (using header)
    blaze_file = os.path.join(params['DRS_CALIB_DB'], hdr['CDBBLAZE'])
    blaze = fits.getdata(blaze_file)
    # normalize blaze
    for order_num in range(blaze.shape[0]):
        blaze[order_num] = blaze[order_num] / np.nanpercentile(blaze[order_num], 95)
    # -------------------------------------------------------------------------
    # get wave solution (using header)
    wave = fits2wave(sp1, hdr)
    # -------------------------------------------------------------------------
    # calculate magic grid for spectrum
    wmin = np.median(wave[first_order])
    wmax = np.median(wave[last_order])

    magic = get_magic_grid(wmin, wmax, dv_grid=1000)

    w = np.zeros_like(magic)
    sp = np.zeros_like(magic)
    middle_nans = np.zeros_like(magic).astype(bool)
    # -------------------------------------------------------------------------
    # plot setup
    plt.close()
    fig, frames = plt.subplots(nrows=3, ncols=1, sharex='all',
                               figsize=[12, 8])
    # -------------------------------------------------------------------------
    # loop around orders in range
    colors = ['blue', 'red', 'cyan']

    min_order = np.zeros(wave.shape[0])
    max_order = np.zeros(wave.shape[0])

    for it, order_num in enumerate(range(first_order, last_order + 1)):
        color = colors[it]
        valid = np.isfinite(sp1[order_num] * blaze[order_num])

        order_name = 79 - order_num

        spl1 = ius(wave[order_num][valid], sp1[order_num][valid], k=1, ext=3)
        spl2 = ius(wave[order_num][valid], blaze[order_num][valid], k=1, ext=3)
        mask = ius(wave[order_num], valid, k=1, ext=3)

        # get nans not on the edge
        wstart = np.min(wave[order_num][valid])
        wend = np.max(wave[order_num][valid])

        wmask = (wave[order_num] > wstart) & (wave[order_num] < wend)

        middlenan = np.isnan(sp1[order_num] * blaze[order_num]) & wmask

        is_middle_nan = ius(wave[order_num], middlenan, k=1, ext=3)

        sp += spl1(magic) * mask(magic)
        w += spl2(magic) * mask(magic)

        middle_nans |= np.array(is_middle_nan(magic)).astype(bool)

        frames[0].plot(wave[order_num], sp1[order_num], color=color,
                       label=f'Order #{order_name}', alpha=0.5, lw=0.5)
        frames[1].plot(wave[order_num], blaze[order_num], color=color,
                       label=f'Order #{order_name}', alpha=0.5, lw=0.5)

        min_order[order_num] = wstart
        max_order[order_num] = wend

    # work out where there is overlap
    omask = np.zeros_like(magic, dtype=bool)
    for order_num in range(wave.shape[0] - 1):
        if min_order[order_num + 1] == 0 or max_order[order_num] == 0:
            continue
        omask |= ((magic < max_order[order_num]) & (magic > min_order[order_num + 1]))

    sp_overlap = sp / w
    sp_no_overlap = sp / w

    sp_overlap[~omask] = np.nan
    sp_no_overlap[omask] = np.nan

    # apply all valid criteria
    sp[middle_nans] = np.nan
    w[middle_nans] = np.nan

    # plot
    frames[0].plot(magic, sp, color='grey', alpha=0.5, label='Combined flux',
                   lw=0.5)
    frames[1].plot(magic, w, color='grey', alpha=0.5, label='Combined flux',
                   lw=0.5)
    frames[2].plot(magic, sp_no_overlap, color='purple', label='no overlap',
                   lw=0.5)
    frames[2].plot(magic, sp_overlap, color='orange', label='overlap', lw=0.5)
    frames[1].set(xlim=[wmin, wmax])

    # -------------------------------------------------------------------------
    # plot labels
    frames[0].set(ylabel='flux [e-]')
    frames[1].set(ylabel='Normalized blaze')
    frames[2].set(ylabel='ratio = spectrum / weight', xlabel='Wavelength [nm]')

    frames[0].legend(loc=0)
    frames[1].legend(loc=0)
    frames[2].legend(loc=0)

    plt.subplots_adjust(hspace=0, wspace=0, left=0.075, right=0.98, top=0.95,
                        bottom=0.075)
    plt.suptitle(hdr['OBJECT'])

    outfile = os.path.join(PLOT_PATH, 's1d.pdf')
    print('Saving to file: ' + outfile)
    plt.savefig(outfile)
    print('Showing graph')
    plt.show(block=True)
    plt.close()


# TCORR
def plot_tcorr_plot(params):
    odocode = '2510303'

    orders = [20, 31]

    limits_a = [[1310, 1320], [1600, 1610]]
    limits_b = [[1314.5, 1315.5], [1600, 1601]]
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
    plt.figure(figsize=(10, 10))
    shape = (2, 5)
    frame1 = plt.subplot2grid(shape, (0, 0), colspan=3)
    frame2 = plt.subplot2grid(shape, (1, 0), colspan=3)
    frame3 = plt.subplot2grid(shape, (0, 3), colspan=2)
    frame4 = plt.subplot2grid(shape, (1, 3), colspan=2)

    frames_right = [frame1, frame2]
    frames_left = [frame3, frame4]

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
                       lw=0.5, label='Extracted', zorder=2)
            # plot telluric
            frame.plot(ordwave, ordtcorr / np.nanmedian(ordtcorr), color='r',
                       lw=0.5, label='Telluric corrected', zorder=3)
            # plot recon
            frame.plot(ordwave, ordrecon / np.nanmedian(ordrecon), color='b',
                       lw=0.5, label='Reconstructed absorption', zorder=4,
                       alpha=0.75)

            frame.plot([0, 0], [1000, 2000],
                       color='orange', lw=0.5, label='Sky model (OH lines)',
                       zorder=1)

            ins = frame.inset_axes([0, 1, 1, 0.1])
            # plot sky
            ins.plot(ordwave, ordskymodel / np.nanmedian(orde2ds),
                     color='orange', lw=0.5, label='Sky model (OH lines)',
                     zorder=1)
            ins.set(xlim=[wavemin, wavemax])
            ins.set_xticklabels([])
            ins.set_yticklabels([])

            # set labels
            if part == 0:
                frame.set(xlabel='Normalized flux')
            frame.set(ylim=[0.4, 1.15], xlim=[wavemin, wavemax])
            frame.set(xlabel='Wavelength [nm]')

            # frame.set_title('Order {0}'.format(order_num), y=1.0, pad=-15)

            if it == 0 and part == 0:
                frame.legend(loc='upper center', bbox_to_anchor=(0.85, 1.2),
                               ncol=5)

            # turn off stuff for zoom
            if part == 1:
                waverange = wavemax - wavemin
                point1 = wavemin + 0.25 * waverange
                point2 = wavemax - 0.25 * waverange
                frame.set_xticks([point1, point2])
                frame.set_xticklabels([str(point1), str(point2)])
                frame.set_yticklabels([])

    plt.subplots_adjust(left=0.075, right=0.98, bottom=0.0725, top=0.9,
                        hspace=0.25, wspace=0.1)
    # -------------------------------------------------------------------------
    # save plot
    outfile = os.path.join(PLOT_PATH, 'tcorr.pdf')
    print('Saving to file: ' + outfile)
    plt.savefig(outfile)
    print('Showing graph')
    plt.show(block=True)
    plt.close()


# TELLU_COV
def plot_tellu_cov_plot(params):
    # get telluric database
    telludbm = drs_database.TelluricDatabase(params)
    telludbm.load_db()

    tellu_table = telludbm.get_tellu_entry('*', key='TELLU_TRANS')

    # group by object name
    uobjnames = np.unique(tellu_table['OBJECT'])
    tau_water = tellu_table['TAU_WATER']
    tau_others = tellu_table['TAU_OTHERS']

    # -------------------------------------------------------------------------
    # plot setup
    plt.close()
    fig, frame = plt.subplots(ncols=1, nrows=1, figsize=(12, 8))

    colors = ['red', 'orange', 'green', 'blue', 'purple', 'black']
    markers = ['o', '+', 's', 'd', '^', 'v']

    fmt = list(itertools.product(markers, colors))

    for it, uobjname in enumerate(uobjnames):

        marker, color = fmt[it]

        if marker in ['+']:
            facecolor = color
            edgecolor = None
        else:
            facecolor = 'None'
            edgecolor = color

        objmask = tellu_table['OBJECT'] == uobjname
        length = np.sum(objmask)

        label = f'{uobjname} N={length}'

        frame.scatter(tau_water[objmask], tau_others[objmask], alpha=0.75,
                      label=label, marker=marker, edgecolor=edgecolor,
                      facecolor=facecolor)

    frame.set(xlabel=r'$\tau$[water]$\approx$water absorption',
              ylabel=r'$ \tau$[dry]$\approx$dry absorption (airmass)')
    frame.legend(title='Hot stars', loc=6, bbox_to_anchor=(1.025, 0.5))
    plt.subplots_adjust(left=0.075, right=0.8, bottom=0.0725, top=0.9,
                        hspace=0.2)
    # -------------------------------------------------------------------------
    # save plot
    outfile = os.path.join(PLOT_PATH, 'tellu_cov.pdf')
    print('Saving to file: ' + outfile)
    plt.savefig(outfile)
    print('Showing graph')
    plt.show(block=True)
    plt.close()


def plot_tellu_airmass_plot(params):
    # get telluric database
    telludbm = drs_database.TelluricDatabase(params)
    telludbm.load_db()

    tellu_table = telludbm.get_tellu_entry('*', key='TELLU_TRANS')

    # group by object name
    uobjnames = np.unique(tellu_table['OBJECT'])
    airmass = tellu_table['AIRMASS']
    tau_others = tellu_table['TAU_OTHERS']

    # -------------------------------------------------------------------------
    # plot setup
    plt.close()
    fig = plt.figure(figsize=(10, 8))
    frame1 = plt.subplot2grid((6, 6), (0, 0), rowspan=5, colspan=6)
    frame2 = plt.subplot2grid((6, 6), (5, 0), rowspan=1, colspan=6)

    colors = ['red', 'orange', 'green', 'blue', 'purple', 'black']
    markers = ['o', '+', 's', 'd', '^', 'v']

    fmt = list(itertools.product(markers, colors))

    for it, uobjname in enumerate(uobjnames):

        marker, color = fmt[it]

        if marker in ['+']:
            facecolor = color
            edgecolor = None
        else:
            facecolor = 'None'
            edgecolor = color

        objmask = tellu_table['OBJECT'] == uobjname
        length = np.sum(objmask)

        label = f'{uobjname} N={length}'

        frame1.scatter(airmass[objmask], tau_others[objmask], alpha=0.75,
                      label=label, marker=marker, edgecolor=edgecolor,
                      facecolor=facecolor)

        frame2.scatter(airmass[objmask], tau_others[objmask] - airmass[objmask],
                       ls='None', marker=marker, edgecolor=edgecolor,
                       facecolor=facecolor)

    limits = [0.9, 2.5]

    frame1.plot(limits, limits, color='k', alpha=0.25)

    frame1.set(xlabel=r'Airmass',
              ylabel=r'$ \tau$[dry]$\approx$dry absorption',
               xlim=limits, ylim=limits)
    frame1.xaxis.tick_top()
    frame1.legend(title='Hot stars', loc=6, bbox_to_anchor=(1.025, 0.4))

    frame2.set(xlabel=r'Airmass', ylabel=r'Difference',
               xlim=limits)

    plt.subplots_adjust(left=0.075, right=0.8, bottom=0.0725, top=0.95,
                        hspace=0.0, wspace=0.0)
    # -------------------------------------------------------------------------
    # save plot
    outfile = os.path.join(PLOT_PATH, 'tellu_airmass.pdf')
    print('Saving to file: ' + outfile)
    plt.savefig(outfile)
    print('Showing graph')
    plt.show(block=True)
    plt.close()


# THERM
def plot_therm_plot(params):
    # file to get
    hashcode = '2510376o'
    # science fibers to correction
    fibers = ['AB']
    # -------------------------------------------------------------------------
    # custom imports only for therm plot
    #   change manually for different instrument
    from apero.base import base
    from apero.core.instruments.spirou import recipe_definitions
    from apero.core.instruments.spirou import file_definitions
    from apero.science.calib import thermal
    from apero.core.constants import ParamDict
    from apero.core.utils import drs_data
    from apero.core.utils import drs_recipe
    # -------------------------------------------------------------------------
    # set the func name
    mainname = __NAME__ + '.plot_therm_plot()'
    # alias pcheck
    pcheck = constants.PCheck(wlog=WLOG)
    # must set a pid (locking required)
    params['PID'] = 'UNKNOWN-PID'
    # must set date now
    params['DATE_NOW'] = base.__now__.fits
    # force no save of products
    if 'INPUTS' not in params:
        params['INPUTS'] = ParamDict()
    params['INPUTS']['NOSAVE'] = True
    # make sure we have a LEAKCORR parameter
    params['INPUTS']['LEAKCORR'] = True
    # set outpath
    params['OUTPATH'] = PLOT_PATH
    params['OBS_DIR'] = NIGHT
    # load pseudo constants
    pconst = constants.pload()
    # get the fiber types needed
    sci_fibers, ref_fiber = pconst.FIBER_KINDS()
    # must process reference fiber first then required fibers
    fibertypes = [ref_fiber] + fibers
    # load the calibration database
    calibdbm = drs_database.CalibrationDatabase(params)
    calibdbm.load_db()
    # get file paths
    pp_file = os.path.join(params['DRS_DATA_WORKING'], NIGHT,
                           '{0}_pp.fits'.format(hashcode))
    # get recipe
    recipe = recipe_definitions.apero_extract
    recipe.plot = drs_recipe.lambda_plot
    # load file
    print('Loading pp image')
    pp_image = fits.getdata(pp_file)
    pp_header = fits.getheader(pp_file)
    # create infile
    infile = file_definitions.pp_file.newcopy(filename=pp_file, params=params,
                                              data=pp_image, header=pp_header)
    # =========================================================================
    # All this comes from apero_extract_spirou.py
    # =========================================================================
    mout1 = mock_extract_part1(params, recipe, infile, fibertypes, calibdbm)
    # get products of mock extract part1
    image, image2, orderps, orderpfiles = mout1[:4]
    orderptimes, props, sprops, bprops = mout1[4:]
    # storage for return / reference fiber usage
    e2dsoutputs = dict()
    eprops_all = dict()
    waves_all = dict()
    # loop around fiber types
    for fiber in fibertypes:

        if fiber == ref_fiber:
            do_thermal = True
        else:
            do_thermal = False

        mout2 = mock_extract_part2(params, recipe, infile, fibertypes,
                                   ref_fiber, fiber, e2dsoutputs, image, image2,
                                   orderps, orderpfiles, orderptimes, props,
                                   sprops, bprops, calibdbm,
                                   do_thermal=do_thermal)
        # get products of mock extract part 2
        eprops, wprops, e2dsoutputs = mout2
        # must save eprops
        eprops_all[fiber] = eprops
        waves_all[fiber] = wprops

    # =========================================================================
    # END OF apero_extract_spirou.py
    # =========================================================================

    # Now we have to manually correct thermal for fibers
    for fiber in fibers:
        # set up
        kwargs, func_name = dict(), mainname
        # get eprops for this fiber
        eprops = eprops_all[fiber]
        blaze = eprops['BLAZE']
        # get the wavemap
        wavemap = waves_all[fiber]['WAVEMAP']
        # get properties from parameter dictionaries / kwargs
        tapas_thres = pcheck(params, 'THERMAL_THRES_TAPAS', 'tapas_thres', kwargs,
                             func_name)
        torder = pcheck(params, 'THERMAL_ORDER', 'torder', kwargs, func_name)
        red_limit = pcheck(params, 'THERMAL_RED_LIMIT', 'red_limit', kwargs,
                           func_name)
        e2ds = pcheck(params, 'E2DS', 'e2ds', kwargs, func_name,
                      paramdict=eprops)
        thermal_file = kwargs.get('thermal_file', None)
        # get the thermal correction file
        tout = thermal.get_thermal(params, pp_header, fiber=fiber,
                                   filename=thermal_file,
                                   kind='THERMALT_E2DS', database=calibdbm,
                                   required=True)
        thermalfile, thermaltime, thermal = tout
        # ----------------------------------------------------------------------
        # load tapas
        tapas, _ = drs_data.load_tapas(params)
        wtapas, ttapas = tapas['wavelength'], tapas['trans_combined']
        # --------------------------------------------------------------------------
        # Method 1: Use tapas (for use with bright targets)
        # --------------------------------------------------------------------------
        # binary mask to be saved; this corresponds to the domain for which
        #    transmission is basically zero and we can safely use the domain
        #    to scale the thermal background. We only do this for wavelength smaller
        #    than "THERMAL_TAPAS_RED_LIMIT" nm as this is the red end of the
        #    TAPAS domain
        # --------------------------------------------------------------------------
        # splining tapas onto the order 49 wavelength grid
        sptapas = mp.iuv_spline(wtapas, ttapas, ext=3, k=1)
        # set torder mask all to False initially
        torder_mask = np.zeros_like(wavemap[torder, :], dtype=bool)
        # get the wave mask
        wavemask = wavemap[torder] < red_limit
        # get the tapas data for these wavelengths
        torder_tapas = sptapas(wavemap[torder, wavemask])
        # find those pixels lower than threshold in tapas
        torder_mask[wavemask] = torder_tapas < tapas_thres
        # we find the median scale between the observation and the thermal
        #    background in domains where there is no transmission
        thermal_torder = thermal[torder, torder_mask]
        image_torder = e2ds[torder, torder_mask]
        ratio = mp.nanmedian(thermal_torder / image_torder)
        # -------------------------------------------------------------------------
        # scale thermal by ratio
        thermal = thermal / ratio
        # get properties from params
        startorder = params['THERMAL_PLOT_START_ORDER']
        # correct data for graph
        rblaze = np.ravel(blaze[startorder:])
        rwave = np.ravel(wavemap[startorder:])
        rimage = np.ravel(e2ds[startorder:])
        rthermal = np.ravel(thermal[startorder:])
        swave = wavemap[torder, torder_mask]
        sthermal = thermal[torder][torder_mask]
        sblaze = blaze[torder][torder_mask]
        # -------------------------------------------------------------------------
        # plot setup
        plt.close()
        fig, frame = plt.subplots(ncols=1, nrows=1, figsize=(10, 4))
        # -------------------------------------------------------------------------
        # plot data
        frame.plot(rwave, rimage / rblaze, color='k', label='input spectrum',
                   lw=0.25)
        frame.plot(rwave, rthermal / rblaze, color='r', label='scaled thermal',
                   lw=0.5)
        frame.plot(swave, sthermal / sblaze, markeredgecolor='b',
                   markerfacecolor='None', marker='o', ls='None',
                   label='background sample region', ms=2, markeredgewidth=0.25)

        # set graph properties
        frame.legend(loc=0)
        frame.set(xlabel='Wavelength [nm]', ylabel='Flux/Blaze',
                  xlim=[2300, 2500], ylim=[-2e-5, 0.0016])
        plt.subplots_adjust(left=0.1, right=0.975, bottom=0.15, top=0.95,
                            hspace=0.2)
        # -------------------------------------------------------------------------
        # save plot
        outfile = os.path.join(PLOT_PATH, f'therm_plot_{fiber}.pdf')
        print('Saving to file: ' + outfile)
        plt.savefig(outfile)
        print('Showing graph')
        plt.show(block=True)
        plt.close()


# LEAK
def plot_leak_plot(params):
    # file to get
    hashcode = '2510286a'
    # science fibers to correction
    fibers = ['AB']
    # -------------------------------------------------------------------------
    # custom imports only for therm plot
    #   change manually for different instrument
    from apero.base import base
    from apero.core.instruments.spirou import recipe_definitions
    from apero.core.instruments.spirou import file_definitions
    from apero.core.constants import ParamDict
    from apero.core.utils import drs_recipe
    # -------------------------------------------------------------------------
    # must set a pid (locking required)
    params['PID'] = 'UNKNOWN-PID'
    # must set date now
    params['DATE_NOW'] = base.__now__.fits
    # force no save of products
    if 'INPUTS' not in params:
        params['INPUTS'] = ParamDict()
    params['INPUTS']['NOSAVE'] = True
    # make sure we have a LEAKCORR parameter
    params['INPUTS']['LEAKCORR'] = True
    # set outpath
    params['OUTPATH'] = PLOT_PATH
    params['OBS_DIR'] = NIGHT
    # load pseudo constants
    pconst = constants.pload()
    # get the fiber types needed
    sci_fibers, ref_fiber = pconst.FIBER_KINDS()
    # must process reference fiber first then required fibers
    fibertypes = [ref_fiber] + fibers
    # load the calibration database
    calibdbm = drs_database.CalibrationDatabase(params)
    calibdbm.load_db()
    # get file paths
    pp_file = os.path.join(params['DRS_DATA_WORKING'], NIGHT,
                           '{0}_pp.fits'.format(hashcode))
    # get recipe
    recipe = recipe_definitions.apero_extract
    recipe.plot = drs_recipe.lambda_plot
    # load file
    print('Loading pp image')
    pp_image = fits.getdata(pp_file)
    pp_header = fits.getheader(pp_file)
    # create infile
    infile = file_definitions.pp_file.newcopy(filename=pp_file, params=params,
                                              data=pp_image, header=pp_header)
    # get colour maps
    cmap1 = matplotlib.cm.get_cmap('inferno').copy()
    cmap1.set_bad(color='green')
    # =========================================================================
    # All this comes from apero_extract_spirou.py
    # =========================================================================
    mout1 = mock_extract_part1(params, recipe, infile, fibertypes, calibdbm)
    # get products of mock extract part1
    image, image2, orderps, orderpfiles = mout1[:4]
    orderptimes, props, sprops, bprops = mout1[4:]
    # storage for return / reference fiber usage
    e2dsoutputs = dict()
    eprops_all = dict()
    waves_all = dict()
    # loop around fiber types
    for fiber in fibertypes:

        if fiber == ref_fiber:
            do_thermal = True
        else:
            do_thermal = False

        mout2 = mock_extract_part2(params, recipe, infile, fibertypes,
                                   ref_fiber, fiber, e2dsoutputs, image, image2,
                                   orderps, orderpfiles, orderptimes, props,
                                   sprops, bprops, calibdbm,
                                   do_thermal=do_thermal)
        # get products of mock extract part 2
        eprops, wprops, e2dsoutputs = mout2
        # must save eprops
        eprops_all[fiber] = eprops
        waves_all[fiber] = wprops

    # =========================================================================
    # END OF apero_extract_spirou.py
    # =========================================================================

    for fiber in fibers:

        # get eprops for this fiber
        eprops = eprops_all[fiber]
        uncorr_e2ds = eprops['UNCORR_E2DS']
        e2ds = eprops['E2DS']
        wavemap = waves_all[fiber]['WAVEMAP']

        # ---------------------------------------------------------------------
        # plot setup
        plt.close()
        fig = plt.figure(figsize=(12, 6))
        size = (2, 2)
        frame1 = plt.subplot2grid(size, (0, 0), colspan=1, rowspan=2)
        frame2 = plt.subplot2grid(size, (0, 1), colspan=1, rowspan=1)
        frame3 = plt.subplot2grid(size, (1, 1), colspan=1, rowspan=1)
        # ---------------------------------------------------------------------

        pp_image_zoom = pp_image[1325:1525, 540:620]

        limit1 = [1720, 1770]
        limit2 = [1729, 1732]

        mask1 = (wavemap[order_num] > limit1[0])
        mask1 &= (wavemap[order_num] < limit1[1])
        mask2 = (wavemap[order_num] > limit2[0])
        mask2 &= (wavemap[order_num] < limit2[1])


        order_num = 35
        # plot pp_image contaminated
        im, _ = imshow_norm(pp_image_zoom, frame1, origin='lower', aspect='auto',
                           interval=ZScaleInterval(), stretch=LogStretch(100),
                           cmap=cmap1, interpolation='None', rasterized=True)

        cbar = _add_colorbar(fig, im, frame1, side='bottom', pad=0.1)
        cbar.set_ticks([0.1, 0.2, 0.4, 0.8, 1.0, 2.0, 4.0, 8.0])

        # plot e2ds order 35 contaminated
        frame2.plot(wavemap[order_num][mask1], uncorr_e2ds[order_num][mask1],
                    color='k', label='input spectrum', lw=0.25)
        # plot e2ds order 35 decontamined
        frame2.plot(wavemap[order_num][mask1], e2ds[order_num][mask1],
                    color='r', label='corrected spectrum', lw=0.25)

        # plot e2ds order 35 contaminated
        frame3.plot(wavemap[order_num][mask2], uncorr_e2ds[order_num][mask2],
                    color='k', label='input spectrum', lw=0.25)
        # plot e2ds order 35 decontamined
        frame3.plot(wavemap[order_num][mask2], e2ds[order_num][mask2],
                    color='r', label='corrected spectrum', lw=0.25)

        # set graph properties
        frame1.set_xticklabels([])
        frame1.set_yticklabels([])
        frame2.legend(loc=0)
        frame2.set(xlabel='Wavelength [nm]', ylabel='Flux', xlim=limit1)

        frame2.xaxis.set_label_position('top')
        frame2.xaxis.tick_top()
        frame2.yaxis.set_label_position('right')
        frame2.yaxis.tick_right()
        frame3.legend(loc=0)
        frame3.set(xlabel='Wavelength [nm]', ylabel='Flux', xlim=limit2)
        frame3.yaxis.set_label_position('right')
        frame3.yaxis.tick_right()

        plt.subplots_adjust(left=0.025, right=0.9, bottom=0.1, top=0.9,
                            hspace=0.025, wspace=0.05)
        # ---------------------------------------------------------------------
        # save plot
        outfile = os.path.join(PLOT_PATH, f'leak_plot_{fiber}.pdf')
        print('Saving to file: ' + outfile)
        plt.savefig(outfile)
        print('Showing graph')
        plt.show(block=True)
        plt.close()


# =============================================================================
# worker functions (private)
# =============================================================================
def _norm_image(image, frame, cmap):
    im, norm = imshow_norm(image, frame, origin='lower', aspect='auto',
                           interval=ZScaleInterval(), stretch=LinearStretch(),
                           cmap=cmap, interpolation='None', rasterized=True)
    return im, norm


def _add_colorbar(fig, im, frame, side='bottom', pad=0.0):
    if side in ['right', 'left']:
        orientation = 'vertical'
    else:
        orientation = 'horizontal'

    divider = make_axes_locatable(frame)
    cax = divider.append_axes(side, '5%', pad=pad)
    cbar = fig.colorbar(im, cax=cax, orientation=orientation)
    cbar.ax.tick_params(labelsize=8)

    return cbar


def poly_region(mask: np.ndarray):
    poly = []

    # get x limits
    xmin, xmax = -0.5, mask.shape[1]

    flatindex = np.where(mask)

    # find the first row
    y_first = flatindex[0].ravel()[0] - 0.5
    # find the first col on the first row
    x_first = flatindex[1].ravel()[0] - 0.5
    # find the last row
    y_last = flatindex[0].ravel()[-1] + 0.5
    # find the low col on the last row
    x_last = flatindex[1].ravel()[-1] + 0.5

    # define shape of poly
    poly.append([x_first, y_first])

    poly.append([xmax, y_first])

    poly.append([xmax, y_last - 1])

    poly.append([x_last, y_last - 1])

    poly.append([x_last, y_last])

    poly.append([xmin, y_last])

    poly.append([xmin, y_first + 1])

    poly.append([x_first, y_first + 1])

    poly.append([x_first, y_first])

    # work out the mixpoint (in y)
    midpoint = [np.mean([xmax, xmin]), np.mean([y_first, y_last])]
    # add x (just middle)

    return Polygon(poly, facecolor='None', edgecolor='w'), midpoint


def fits2wave(image, header):
    """
    Get the wave solution from the header using a filename
    """
    # size of the image
    nbypix, nbxpix = image.shape
    # get the keys with the wavelength polynomials
    wave_hdr = header['WAVE0*']
    # concatenate into a numpy array
    wave_poly = np.array([wave_hdr[i] for i in range(len(wave_hdr))])
    # get the number of orders
    nord = header['WAVEORDN']
    # get the per-order wavelength solution
    wave_poly = wave_poly.reshape(nord, len(wave_poly) // nord)
    # project polynomial coefficiels
    wavesol = np.zeros_like(image)
    # loop around orders
    for order_num in range(nord):
        ordwave = np.polyval(wave_poly[order_num][::-1], np.arange(nbxpix))
        wavesol[order_num] = ordwave
    # return wave grid
    return wavesol


def get_magic_grid(wave0: Union[np.ndarray, int] = 1500,
                   wave1: Union[np.ndarray, int] = 1800,
                   dv_grid=0.5):
    # default for the function is 500 m/s
    # the arithmetic is a but confusing here, you first find how many
    # elements you have on your grid, then pass it to an exponential
    # the first element is exactely wave0, the last element is NOT
    # exactly wave1, but is very close and is set to get your exact
    # step in velocity
    len_magic = int(np.ceil(np.log(wave1 / wave0) * np.array(cc.c.value) / dv_grid))
    magic_grid = np.exp(np.arange(len_magic) / len_magic * np.log(wave1 / wave0)) * wave0
    return magic_grid


# =============================================================================
# extraction functions
# =============================================================================
def mock_extract_part1(params, recipe, infile, fibertypes, calibdbm):
    # -------------------------------------------------------------------------
    # custom imports only for therm plot
    #   change manually for different instrument
    from apero.science.calib import gen_calib, shape
    from apero.science import extract
    # get header from infile
    header = infile.header
    # =========================================================================
    # START OF apero_extract_spirou.py
    # =========================================================================
    # -------------------------------------------------------------------------
    # Load shape components
    sprops = shape.get_shape_calibs(params, header, database=calibdbm)
    # -------------------------------------------------------------------------
    # Correction of file
    props, image = gen_calib.calibrate_ppfile(params, recipe, infile,
                                              database=calibdbm)
    # -------------------------------------------------------------------------
    # Load and straighten order profiles
    sargs = [infile, fibertypes, sprops]
    oout = extract.order_profiles(params, recipe, *sargs, database=calibdbm)
    orderps, orderpfiles, orderptimes = oout
    # -------------------------------------------------------------------------
    # Apply shape transformations
    # log progress (straightening orderp)
    WLOG(params, 'info', textentry('40-016-00004'))
    # straighten image
    image2 = shape.ea_transform(params, image, sprops['SHAPEL'],
                                dxmap=sprops['SHAPEX'], dymap=sprops['SHAPEY'])
    # -------------------------------------------------------------------------
    # Calculate Barycentric correction
    bprops = extract.get_berv(params, infile, header)

    # return everything needed for mock_extract_part2 (inside fiber loop)
    return [image, image2, orderps, orderpfiles, orderptimes, props, sprops,
            bprops]


def mock_extract_part2(params, recipe, infile, fibertypes, ref_fiber,
                       fiber, e2dsoutputs, image, image2, orderps, orderpfiles,
                       orderptimes, props, sprops, bprops, calibdbm,
                       do_thermal=True, do_leak=True):
    # -------------------------------------------------------------------------
    # custom imports only for therm plot
    #   change manually for different instrument
    from apero.science.calib import shape, wave, localisation
    from apero.science.calib import flat_blaze, leak, thermal
    from apero.science import extract
    from apero.core.constants import ParamDict
    # get header from infile
    header = infile.header
    # set function name
    mainname = 'mock_extract_part2'
    # =========================================================================
    # FROM apero_extract_spirou.py
    # =========================================================================
    # log process: processing fiber
    wargs = [fiber, ', '.join(fibertypes)]
    WLOG(params, 'info', textentry('40-016-00014', args=wargs))
    # ------------------------------------------------------------------
    # get reference fiber data
    ref_key = 'E2DS_{0}'.format(ref_fiber)
    # if we have reference data populate ref_e2ds
    if ref_key in e2dsoutputs:
        ref_e2ds = e2dsoutputs[ref_key].data
    # otherwise this is set to None - and we cannot use it
    else:
        ref_e2ds = None
    # ------------------------------------------------------------------
    # check forcing reference wave solution
    mwave = False
    if 'FORCE_REF_WAVE' in params['INPUTS']:
        mwave = params['INPUTS']['FORCE_REF_WAVE']
    # get the wave solution
    wprops = wave.get_wavesolution(params, recipe, header,
                                   fiber=fiber, ref=mwave,
                                   database=calibdbm,
                                   nbpix=image.shape[1])
    # --------------------------------------------------------------
    # load the localisation properties for this fiber
    lprops = localisation.get_coefficients(params, header, fiber=fiber,
                                           merge=True,
                                           database=calibdbm)
    # get the localisation center coefficients for this fiber
    lcoeffs = lprops['CENT_COEFFS']
    # shift the coefficients
    lcoeffs2 = shape.ea_transform_coeff(image2, lcoeffs,
                                        sprops['SHAPEL'])
    # --------------------------------------------------------------
    # load the flat file for this fiber
    fout = flat_blaze.get_flat(params, header, fiber, database=calibdbm)
    # --------------------------------------------------------------
    # load the blaze file for this fiber
    bout = flat_blaze.get_blaze(params, header, fiber,
                                database=calibdbm)
    # add blaze and flat to parameter dictionary
    fbprops = ParamDict()
    fbprops['FLAT'] = fout[2]
    fbprops['FLATFILE'] = fout[0]
    fbprops['FLATTIME'] = fout[1]
    fbprops['BLAZE'] = bout[2]
    fbprops['BLAZEFILE'] = bout[0]
    fbprops['BLAZETIME'] = bout[1]
    # add keys
    keys = ['FLAT', 'FLATFILE', 'FLATTIME', 'BLAZE', 'BLAZEFILE',
            'BLAZETIME']
    fbprops.set_sources(keys, mainname)
    # get the number of frames used
    nframes = 1
    # --------------------------------------------------------------
    # get the order profile for this fiber
    lprops['ORDERP'] = orderps[fiber]
    lprops['ORDERPFILE'] = orderpfiles[fiber]
    lprops['ORDERPTIME'] = orderptimes[fiber]
    lprops.set_sources(['ORDERP', 'ORDERPFILE', 'ORDERPTIME'], mainname)
    # --------------------------------------------------------------
    # log progress: extracting image
    WLOG(params, 'info', textentry('40-016-00011'))
    # extract spectrum
    eprops = extract.extract2d(params, image2, lprops['ORDERP'],
                               lcoeffs2, nframes, props, fiber=fiber)
    # leak correction
    if do_leak:
        eprops = leak.manage_leak_correction(params, recipe, eprops,
                                             infile, fiber, ref_e2ds)
    # flat correction for e2dsff
    eprops = extract.flat_blaze_correction(eprops, fbprops['FLAT'],
                                           fbprops['BLAZE'])
    # do not do thermal correction here other than for reference fiber
    #    (we will do it afterwards)
    if do_thermal:
        eprops = thermal.thermal_correction(params, recipe, header,
                                            props, eprops, fiber=fiber,
                                            database=calibdbm)
    # --------------------------------------------------------------
    s1dextfile = params['EXT_S1D_INTYPE']
    # create 1d spectra (s1d) of the e2ds file
    sargs = [wprops['WAVEMAP'], eprops[s1dextfile], eprops['BLAZE']]
    swprops = extract.e2ds_to_s1d(params, recipe, *sargs,
                                  wgrid='wave', fiber=fiber,
                                  s1dkind=s1dextfile)
    svprops = extract.e2ds_to_s1d(params, recipe, *sargs,
                                  wgrid='velocity', fiber=fiber,
                                  s1dkind=s1dextfile)

    # --------------------------------------------------------------
    # Quality control
    qc_params, passed = extract.qc_extraction(params, eprops)
    # --------------------------------------------------------------
    # add files to outputs
    # --------------------------------------------------------------
    fargs = [params, recipe, infile, [infile.basename], False, fiber,
             props, lprops, wprops, eprops, bprops,
             swprops, svprops, sprops, fbprops, qc_params]
    outfiles = extract.write_extraction_files(*fargs)
    e2dsfile, e2dsfffile = outfiles

    ekeys = ['E2DS', 'E2DSFF']
    efiles = [e2dsfile, e2dsfffile]
    # loop around keys to add
    for key, efile in zip(ekeys, efiles):
        # construct output key
        outkey = '{0}_{1}'.format(key, fiber)
        # copy file to dictionary
        e2dsoutputs[outkey] = efile.completecopy(efile)

    # return product for use afterwards
    return eprops, wprops, e2dsoutputs


# =============================================================================
# backmap functions
# =============================================================================
def _do_bad_pix(params, flat_image, dark_image):
    # ------------------------------------------------------------------
    # Normalise flat and median of flat
    # ------------------------------------------------------------------
    flat_med, flat_image = badpix.normalise_median_flat(params, flat_image)
    # ------------------------------------------------------------------
    # Locate bad pixels
    # ------------------------------------------------------------------
    # Locate bad pixels from dark and flat
    bargs = [flat_image, flat_med, dark_image]
    badpixelmap_a, bstats_a = badpix.locate_bad_pixels(params, *bargs)
    # Locate bad pixels from full detector flat
    bargs = [flat_image]
    badpixelmap_b, bstats_b = badpix.locate_bad_pixels_full(params, *bargs)

    # ------------------------------------------------------------------
    # Combine bad pixel masks
    # ------------------------------------------------------------------
    bad_pixel_map = badpixelmap_a | badpixelmap_b
    # total number of bad pixels
    btotal = (mp.nansum(bad_pixel_map) / bad_pixel_map.size) * 100
    # log result
    WLOG(params, '', textentry('40-012-00007', args=[btotal]))

    # ------------------------------------------------------------------
    # Flip images
    # ------------------------------------------------------------------
    if params['INPUT_FLIP_IMAGE']:
        # flip flat
        flat_image1 = drs_image.flip_image(params, flat_image)
        # flip bad pixel map
        bad_pixel_map1 = drs_image.flip_image(params, bad_pixel_map)
    else:
        flat_image1, bad_pixel_map1 = flat_image, bad_pixel_map

    # ------------------------------------------------------------------
    # Resize image
    # ------------------------------------------------------------------
    if params['INPUT_RESIZE_IMAGE']:
        # get resize size
        sargs = dict(xlow=params['IMAGE_X_LOW'],
                     xhigh=params['IMAGE_X_HIGH'],
                     ylow=params['IMAGE_Y_LOW'],
                     yhigh=params['IMAGE_Y_HIGH'])
        # resize flat
        flat_image1 = drs_image.resize(params, flat_image1, **sargs)
        # resize bad pixel map
        bad_pixel_map1 = drs_image.resize(params, bad_pixel_map1, **sargs)
    else:
        flat_image1 = np.array(flat_image)
        bad_pixel_map1 = np.array(bad_pixel_map)

    backmask, backest = _create_background_map(params, flat_image1,
                                               bad_pixel_map1)
    return flat_image1, backmask, backest


def _create_background_map(params, image, badpixmask):
    # get constants
    width = params['BKGR_BOXSIZE']
    percent = params['BKGR_PERCENTAGE']
    csize = params['BKGR_MASK_CONVOLVE_SIZE']
    nbad = params['BKGR_N_BAD_NEIGHBOURS']
    # set image bad pixels to NaN
    image0 = np.array(image)
    badmask = np.array(badpixmask, dtype=bool)
    image0[badmask] = np.nan
    # image that will contain the background estimate
    backest = np.zeros_like(image0)
    # we slice the image in ribbons of width "width".
    # The slicing is done in the cross-dispersion direction, so we
    # can simply take a median along the "fast" dispersion to find the
    # order profile. We pick width to be small enough for the orders not
    # to show a significant curvature within w pixels
    for x_it in range(0, image0.shape[1], width):
        # ribbon to find the order profile
        ribbon = mp.nanmedian(image0[:, x_it:x_it + width], axis=1)

        for y_it in range(image0.shape[0]):
            # we perform a running Nth percentile filter along the
            # order profile. The box of the filter is w. Note that it could
            # differ in princile from w, its just the same for the sake of
            # simplicity.
            ystart = y_it - width // 2
            yend = y_it + width // 2
            if ystart < 0:
                ystart = 0
            if yend > image0.shape[0] - 1:
                yend = image0.shape[0] - 1
            # background estimate
            backest_pix = np.nanpercentile(ribbon[ystart:yend], percent)
            backest[y_it, x_it: x_it + width] = backest_pix
    # the mask is the area that is below then Nth percentile threshold
    with warnings.catch_warnings(record=True) as _:
        backmask = np.array(image0 < backest, dtype=float)
    # we take advantage of the order geometry and for that the "dark"
    # region of the array be continuous in the "fast" dispersion axis
    nribbon = convolve2d(backmask, np.ones([1, csize]), mode='same')
    # we remove from the binary mask all isolated (within 1x7 ribbon)
    # dark pixels
    backmask[nribbon == 1] = 0
    # If a pixel has 3 or more "dark" neighbours, we consider it dark
    # regardless of its initial value
    backmask[nribbon >= nbad] = 1
    # return the background mask
    return backmask, backest


# =============================================================================
# Start of code
# =============================================================================
if __name__ == '__main__':

    _params = constants.load()

    if 'SIZE_GRID' in PLOTS:
        plot_size_grid(_params)
    if 'FIBER_LAYOUT' in PLOTS:
        plot_fiber_layout(_params)
    if 'RAW_FEATURES' in PLOTS:
        plot_raw_features(_params)
    if 'PP_FEATURES' in PLOTS:
        plot_pp_features(_params)
    if 'PP_CARTOON' in PLOTS:
        plot_pp_cartoon()
    if 'BADMAP' in PLOTS:
        plot_badpix_plot(_params)
    if 'BACKMAP' in PLOTS:
        plot_backmap_plot(_params)
    if 'FLATBLAZE' in PLOTS:
        plot_flatblaze_plot(_params)
    if 'E2DS' in PLOTS:
        plot_e2ds_plot(_params)
    if 'S1D' in PLOTS:
        plot_s1d_plot(_params)
    if 'TCORR' in PLOTS:
        plot_tcorr_plot(_params)
    if 'TELLU_COV' in PLOTS:
        plot_tellu_cov_plot(_params)
    if 'AIRMASS' in PLOTS:
        plot_tellu_airmass_plot(_params)
    if 'THERM' in PLOTS:
        plot_therm_plot(_params)
    if 'LEAK' in PLOTS:
        plot_leak_plot(_params)

# =============================================================================
# End of code
# =============================================================================
