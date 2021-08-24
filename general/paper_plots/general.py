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

from apero import lang
from apero.core import constants
from apero.core.core import drs_log
from apero.core.core import drs_database
from apero.science import preprocessing as prep
from apero.io import drs_image
from apero.core.instruments.spirou import file_definitions


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
# PLOT_PATH = '/data/spirou/drs-data/misc/paper_plots'
PLOT_PATH = '/scratch2/drs-data/misc/paper_plots'
# define plots and append those we want
PLOTS = []
# PLOTS.append('SIZE_GRID')
PLOTS.append('RAW_FEATURES')
# PLOTS.append('BADMAP')


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

# =============================================================================
# End of code
# =============================================================================
