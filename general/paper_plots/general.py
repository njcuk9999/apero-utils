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
from astropy.visualization import imshow_norm, ZScaleInterval, LinearStretch
import glob
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
# define plots and append those we want
PLOTS = []
# PLOTS.append('SIZE_GRID')
# PLOTS.append('RAW_FEATURES')
# PLOTS.append('BADMAP')
# PLOTS.append('BACKMAP')
PLOTS.append('FLATBLAZE')


# =============================================================================
# PLOT functions
# =============================================================================
# SIZE_GRID
def plot_size_grid(params):

    odocode = '2510288a'
    hashcode = '3444961B5D'
    # get file paths
    raw_file = os.path.join(params['DRS_DATA_RAW'], NIGHT,
                            '{0}.fits'.format(odocode))
    pp_file = os.path.join(params['DRS_DATA_REDUC'], NIGHT,
                           '{0}_pp.fits'.format(hashcode))
    e2ds_file = os.path.join(params['DRS_DATA_REDUC'], NIGHT,
                             '{0}_pp_e2dsff_AB.fits'.format(hashcode))
    # get
    print('Loading raw image')
    raw_image = fits.getdata(raw_file)
    print('Loading pp image')
    pp_image = fits.getdata(pp_file)
    print('Loading E2DS image')
    e2ds_image = fits.getdata(e2ds_file)

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
    fig = plt.figure(figsize=(12, 12))
    size = (2, 2)
    frame1 = plt.subplot2grid(size, (0, 0))
    frame2 = plt.subplot2grid(size, (0, 1))
    frame3 = plt.subplot2grid(size, (1, 0), colspan=2)

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

    plt.subplots_adjust(wspace=0.05, hspace=0.05,
                        left=0.01, right=0.99, top=0.975, bottom=0.025)
    # -------------------------------------------------------------------------
    outfile = os.path.join(PLOT_PATH, 'size_grid.pdf')
    print('Saving to file: ' + outfile)
    plt.savefig(outfile, dpi=500)
    print('Showing graph')
    plt.show()
    plt.close()


# RAW_FEATURES
def plot_raw_features(params):
    # set file to use
    odocode = '2510288a'
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
    odocode = '2510288a'
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

        norm = np.nanmedian(e2ds_image[order_num])

        frame0.plot(e2ds_image[order_num]/norm, label=f'Order={order_num}',
                    color=colors[it], alpha=0.5)
        frame0.plot(blaze_image[order_num]/norm, color=colors[it], linestyle='--')
        frame1.plot(flat_image[order_num],color=colors[it], alpha=0.75)

    frame0.legend(loc=0)
    frame0.set(xlabel='Pixel number', ylabel='Normalized E2DS Flux',
               xlim=[0, 4088])
    frame1.set(xlabel='Pixel number', ylabel='Residual E2DS-Blaze Flux',
               xlim=[0, 4088])
    # adjust edges of figure
    plt.subplots_adjust(wspace=0, hspace=0.2, left=0.075, right=0.99,
                        bottom=0.075, top=0.99)
    # -------------------------------------------------------------------------
    # save plot
    outfile = os.path.join(PLOT_PATH, 'flat_blaze.pdf')
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
    if 'BADMAP' in PLOTS:
        plot_badpix_plot(params)
    if 'BACKMAP' in PLOTS:
        plot_backmap_plot(params)
    if 'FLATBLAZE' in PLOTS:
        plot_flatblaze_plot(params)


# =============================================================================
# End of code
# =============================================================================
