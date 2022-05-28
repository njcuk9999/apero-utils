#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plots for the apero drs paper

Created on 2021-08-01

@author: cook
"""
import matplotlib
matplotlib.use('Qt5Agg')
from astropy.io import fits
from astropy.table import Table
from astropy.visualization import imshow_norm, ZScaleInterval, LinearStretch
import glob
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Polygon
import matplotlib.patheffects as PathEffects
import numpy as np
import os
from scipy.signal import convolve2d
import warnings

from apero import lang
from apero.core import constants
from apero.core.core import drs_log
from apero.core.core import drs_database
from apero.science import preprocessing as prep
from apero.io import drs_image
from apero.core.instruments.spirou import file_definitions
from apero.science.calib import badpix
from apero.core import math as mp
from apero.recipes.spirou import apero_extract_spirou

# =============================================================================
# Define variables
# =============================================================================
# Get Logging function
WLOG = drs_log.wlog
# Get the text types
textentry = lang.textentry
# Raw prefix
RAW_PREFIX = file_definitions.raw_prefix
# get the object database
ObjectDatabase = drs_database.ObjectDatabase
# define the night of data we want to use
NIGHT = '2020-08-31'
# define where we want to save plots
PLOT_PATH = '/scratch2/spirou/drs-data/misc/paper_plots'
# define Y, J, H, K
BANDS = dict()
# BANDS['Y'] = [944.126, 1108.771]            # UKIRT Y
# BANDS['J'] = [1080.647, 1406.797]           # 2MASS J
# BANDS['H'] = [1478.738, 1823.102]           # 2MSSS H
# BANDS['K'] = [1954.369, 2344.240]           # 2MASS Ks
BANDS['$Y$'] = [9386.00/10, 11134.00/10]
BANDS['$J$'] = [11535.86/10, 13544.22/10]     # MKO
BANDS['$H$'] = [14628.97/10, 18085.44/10]     # MKO
BANDS['$K_{s}$'] = [19577.92/10, 23431.05/10]     # MKO



# define plots and append those we want
PLOTS = []
# PLOTS.append('SIZE_GRID')
# PLOTS.append('RAW_FEATURES')
# PLOTS.append('PP_FEATURES')
PLOTS.append('BADMAP')
# PLOTS.append('BACKMAP')
# PLOTS.append('FLATBLAZE')
# PLOTS.append('E2DS')
PLOTS.append('TCORR')

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
    fig = plt.figure(figsize=(12, 14))
    size = (5, 4)
    frame1 = plt.subplot2grid(size, (0, 0), colspan=2, rowspan=2)
    frame2 = plt.subplot2grid(size, (0, 2), colspan=2, rowspan=2)
    frame3 = plt.subplot2grid(size, (2, 0), colspan=4, rowspan=2)
    frame4 = plt.subplot2grid(size, (4, 0), colspan=4, rowspan=1)

    cmap = matplotlib.cm.get_cmap('inferno').copy()
    cmap.set_bad(color='green')
    # -------------------------------------------------------------------------
    # top left raw image
    im, norm = _norm_image(raw_image, frame1, cmap)
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
    im, norm = _norm_image(image3, frame2, cmap)
    # add labels
    frame2.set(xlim=(-4, 4092), ylim=(0-250, 4096-250))
    frame2.tick_params(axis='both', which='both', bottom=False, top=False,
                       left=False, right=False, labelleft=False,
                       labelbottom=False)
    frame2.set_title('pre-processed, flipped, resized (3100x4088)', loc='left',
                     x=0.05, y=0.95, pad=-14,
                     color='black', backgroundcolor='white')
    # -------------------------------------------------------------------------
    # bottom: e2ds
    im, norm = _norm_image(e2ds_image, frame3, cmap)
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

        txt.set_path_effects([PathEffects.withStroke(linewidth=1,
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
        frame4.fill_betweenx([-0.2*maxflux, 1.05*maxflux], band[0], band[1],
                             color='0.5', alpha=0.5, zorder=0)
        txt = frame4.text(np.mean(band), 0.75*maxflux, bandname,
                          color='w', zorder=10, fontsize=16, ha='center')
        txt.set_path_effects([PathEffects.withStroke(linewidth=1,
                                                     foreground='k')])
    # plot cosmetics
    frame4.axes.yaxis.set_ticks([])
    frame4.set(xlim=[950, 2500], ylim=[0, 1.05*maxflux], xlabel='Wavelength [nm]')
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
    # storage
    snr_hotpix, rms_list = [], []
    shiftdx, shiftdy = 0, 0
    # do this iteratively as if there is a shift need to re-workout QC
    for iteration in range(2):
        # get pass condition
        cout = prep.test_for_corrupt_files(params, image, hotpixels)
        snr_hotpix, rms_list = cout[0], cout[1]
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
    panel.text(0.5*(cut1area[1] + cut1area[0]), cut1area[3] + textheight,
               'D', color=rcolor, backgroundcolor='white')
    # -------------------------------------------------------------------------
    # add rectangle holes
    rcut2 = patches.Rectangle((cut2area[0], cut2area[2]),
                              cut2area[1] - cut2area[0],
                              cut2area[3] - cut2area[2],
                              linewidth=1, edgecolor=rcolor, facecolor='none')
    panel.add_patch(rcut2)
    # add text
    panel.text(0.5*(cut2area[1] + cut2area[0]), cut2area[3] + textheight,
               'E', color=rcolor, backgroundcolor='white')
    # -------------------------------------------------------------------------
    # add rectangle zoom 1
    rzoom0 = patches.Rectangle((zoom0area[0], zoom0area[2]),
                               zoom0area[1] - zoom0area[0],
                               zoom0area[3] - zoom0area[2],
                               linewidth=1, edgecolor=rcolor, facecolor='none')
    panel.add_patch(rzoom0)
    # add text
    panel.text(0.5*(zoom0area[1] + zoom0area[0]), zoom0area[3] + textheight,
               'A', color=rcolor, backgroundcolor='white')
    # -------------------------------------------------------------------------
    # add rectangle zoom 1
    rzoom1 = patches.Rectangle((zoom1area[0], zoom1area[2]),
                               zoom1area[1] - zoom1area[0],
                               zoom1area[3] - zoom1area[2],
                               linewidth=1, edgecolor=rcolor, facecolor='none')
    panel.add_patch(rzoom1)
    # add text
    panel.text(0.5*(zoom1area[1] + zoom1area[0]), zoom1area[3] + textheight,
               'B', color=rcolor, backgroundcolor='white')
    # -------------------------------------------------------------------------
    # add rectangle zoom 2
    rzoom2 = patches.Rectangle((zoom2area[0], zoom2area[2]),
                               zoom2area[1] - zoom2area[0],
                               zoom2area[3] - zoom2area[2],
                               linewidth=1, edgecolor=rcolor, facecolor='none')
    panel.add_patch(rzoom2)
    # add text
    panel.text(0.5*(zoom2area[1] + zoom2area[0]), zoom2area[3] + textheight,
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

    im, norm = imshow_norm(dark_image, frame0, origin='lower', aspect='auto',
                           interval=ZScaleInterval(), stretch=LinearStretch(),
                           cmap=cmap1, interpolation='None', rasterized=True)
    im, norm = imshow_norm(bad_image_full, frame1, origin='lower', aspect='auto',
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
             100 * np.sum(bad_image_full)/np.product(bad_image.shape)]
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
    # -------------------------------------------------------------------------
    # get colour maps
    cmap1 = matplotlib.cm.get_cmap('inferno').copy()
    cmap2 = matplotlib.cm.get_cmap('Greys_r').copy()
    cmap3 = matplotlib.cm.get_cmap('Greys').copy()
    # -------------------------------------------------------------------------
    # plot setup
    plt.close()
    fig, frames = plt.subplots(figsize=(18, 6), ncols=3, nrows=1)
    frame0 = frames[0]
    frame1 = frames[1]
    frame2 = frames[2]
    # -------------------------------------------------------------------------
    # three imshow norm plots
    im, norm = imshow_norm(flat_image1, frame0, origin='lower', aspect='auto',
                           interval=ZScaleInterval(), stretch=LinearStretch(),
                           cmap=cmap2, interpolation='None', rasterized=True)
    im, norm = imshow_norm(backest, frame1, origin='lower', aspect='auto',
                           cmap=cmap1, interpolation='None', rasterized=True)
    im, norm = imshow_norm(backmask, frame2, origin='lower', aspect='auto',
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
    fig = plt.figure(figsize=(12, 8))


    frame0 = plt.subplot2grid((3, 1), (0, 0), rowspan=2, colspan=1)
    frame1 = plt.subplot2grid((3, 1), (2, 0), rowspan=1, colspan=1)


    for it, order_num in enumerate(orders):

        order_name = 79 - order_num

        # norm = np.nanmedian(e2ds_image[order_num])

        frame0.plot(e2ds_image[order_num], label=f'Order #{order_name}',
                    color=colors[it], alpha=0.5)
        frame0.plot(blaze_image[order_num], color=colors[it], linestyle='--')
        frame1.plot(flat_image[order_num],color=colors[it], alpha=0.75)

    frame0.legend(loc=0)
    frame0.set(xlabel='Pixel number', ylabel='E2DS Flux',
               xlim=[0, 4088])
    frame1.set(xlabel='Pixel number', ylabel='Residual E2DS-Blaze Flux',
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
    im, norm = imshow_norm(pp_image, frame0, origin='lower', aspect='auto',
                           interval=ZScaleInterval(), stretch=LinearStretch(),
                           cmap=cmap1, interpolation='None', rasterized=True)
    im, norm = imshow_norm(e2dsll_image, frame1, origin='lower', aspect='auto',
                           cmap=cmap1, interpolation='None', rasterized=True)
    im, norm = imshow_norm(e2ds_image, frame2, origin='lower', aspect='auto',
                           cmap=cmap1, interpolation='None', rasterized=True)

    im, norm = imshow_norm(pp_image1, frame3, origin='lower', aspect='auto',
                           interval=ZScaleInterval(), stretch=LinearStretch(),
                           cmap=cmap1, interpolation='None', rasterized=True)
    im, norm = imshow_norm(e2dsll_image1, frame4, origin='lower', aspect='auto',
                           cmap=cmap1, interpolation='None', rasterized=True)
    im, norm = imshow_norm(e2ds_image1, frame5, origin='lower', aspect='auto',
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



# TCORR
def plot_tcorr_plot(params):
    odocode = '2510303'

    orders = [20, 31]

    limits = [[1310, 1320], [1600, 1610]]
    # -------------------------------------------------------------------------
    # get filename
    e2dsfile = os.path.join(params['DRS_DATA_OUT'], NIGHT, '{0}e.fits'.format(odocode))
    tcorrfile = os.path.join(params['DRS_DATA_OUT'], NIGHT, '{0}t.fits'.format(odocode))

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
    fig, frames = plt.subplots(ncols=1, nrows=len(orders), figsize=(12, 12))

    # loop around orders
    for it, order_num in enumerate(orders):
        # get this iterations frame
        frame = frames[it]
        # get order values
        ordwave = wave[order_num]
        orde2ds = (e2dsffab/blaze)[order_num]
        ordtcorr = (tcorrab/blaze)[order_num]
        ordrecon = recon[order_num]
        ordskymodel = skymodel[order_num]

        # get first and last point on e2ds
        imask = np.where(np.isfinite(orde2ds))

        # get wave limits for first and last poitns of e2ds
        if limits[it] is None:
            wavemin = ordwave[int(np.min(imask))]
            wavemax = ordwave[int(np.max(imask))]
        else:
            wavemin = limits[it][0]
            wavemax = limits[it][1]


        # plot e2dsff
        frame.plot(ordwave, orde2ds/np.nanmedian(orde2ds), color='k', lw=0.5,
                   label='Extracted', zorder=2)
        # plot telluric
        frame.plot(ordwave, ordtcorr/np.nanmedian(ordtcorr), color='r', lw=0.5,
                   label='Telluric corrected', zorder=3)
        # plot recon
        frame.plot(ordwave, ordrecon/np.nanmedian(ordrecon), color='b', lw=0.5,
                   label='Reconstructed absorption', zorder=4, alpha=0.75)

        frame.plot([0, 0], [1000, 2000],
                   color='orange', lw=0.5, label='Sky model (OH lines)',
                   zorder=1)

        ins = frame.inset_axes([0, 1, 1, 0.1])
        # plot sky
        ins.plot(ordwave, ordskymodel/np.nanmedian(orde2ds),
                 color='orange', lw=0.5, label='Sky model (OH lines)',
                 zorder=1)
        ins.set(xlim=[wavemin, wavemax])
        ins.set_xticklabels([])
        ins.set_yticklabels([])

        # set labels
        frame.set(ylim=[0.4, 1.15], ylabel='Normalized flux',
                  xlim=[wavemin, wavemax])

        frame.set(xlabel='Wavelength [nm]')

        # frame.set_title('Order {0}'.format(order_num), y=1.0, pad=-15)

        if it == 0:
            frame.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=5)

    plt.subplots_adjust(left=0.075, right=0.98, bottom=0.0725, top=0.9,
                        hspace=0.2)
    # -------------------------------------------------------------------------
    # save plot
    outfile = os.path.join(PLOT_PATH, 'tcorr.pdf')
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


def _add_colorbar(fig, im, frame, side='bottom'):

    if side in ['right', 'left']:
        orientation = 'vertical'
    else:
        orientation = 'horizontal'

    divider = make_axes_locatable(frame)
    cax = divider.append_axes(side, '5%', pad=0)
    cbar = fig.colorbar(im, cax=cax, orientation=orientation)
    cbar.ax.tick_params(labelsize=8)


def poly_region(mask: np.ndarray):
    poly = []

    # get x limits
    xmin, xmax = -0.5, mask.shape[1]

    flatindex = np.where(mask)
    index = np.indices(mask.shape)

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
    midpoint = []
    # add x (just middle)
    midpoint.append(np.mean([xmax, xmin]))
    midpoint.append(np.mean([y_first, y_last]))

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
        ordwave = np.polyval(wave_poly[order_num][::-1],np.arange(nbxpix))
        wavesol[order_num] = ordwave
    # return wave grid
    return wavesol


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

    params = constants.load()

    if 'SIZE_GRID' in PLOTS:
        plot_size_grid(params)
    if 'RAW_FEATURES' in PLOTS:
        plot_raw_features(params)
    if 'PP_FEATURES' in PLOTS:
        plot_pp_features(params)
    if 'BADMAP' in PLOTS:
        plot_badpix_plot(params)
    if 'BACKMAP' in PLOTS:
        plot_backmap_plot(params)
    if 'FLATBLAZE' in PLOTS:
        plot_flatblaze_plot(params)
    if 'E2DS' in PLOTS:
        plot_e2ds_plot(params)
    if 'TCORR' in PLOTS:
        plot_tcorr_plot(params)

# =============================================================================
# End of code
# =============================================================================
