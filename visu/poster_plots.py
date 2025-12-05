#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-05-12 at 09:41

@author: cook
"""
import os
from typing import Union

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.visualization import ZScaleInterval
import matplotlib.dates as mdates
from matplotlib import colormaps
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as mticker

# =============================================================================
# Define variables
# =============================================================================
RAW_IMAGE = '/scratch2/spirou/drs-data/common/minidata2/2020-08-31/2510301o.fits'

S1D = '/scratch2/spirou/drs-data/spirou_minidata2_07286_jupiter/out/2020-08-31/2510301s.fits'

RDB = '/scratch2/spirou/misc/poster_CFHTUM25/lbl_GL699_GL699_spirou_offline.rdb'


# -----------------------------------------------------------------------------
BANDS = [[1082, 1085], [1600, 1604], [2164, 2169]]

PLOTS = ['s1d']

FONTSIZE = 16

INFERNO = colormaps['inferno']

PLOT_DIR = '/scratch2/spirou/misc/poster_CFHTUM25'

# =============================================================================
# Define functions
# =============================================================================
def raw_image():
    interval = ZScaleInterval()  # For DS9-like autoscaling
    # interval = MinMaxInterval()  # For simple min-max scaling
    # load image
    data = fits.getdata(RAW_IMAGE)
    data = np.rot90(data)
    # get limits
    vmin, vmax = interval.get_limits(data)
    # Set figure background to black
    plt.figure(figsize=(32, 32), facecolor='black')
    # plot image
    plt.imshow(data, aspect='auto', origin='lower', vmin=vmin, vmax=vmax,
               cmap='inferno')
    # Turn off axes
    plt.axis('off')
    # Set axes background to black
    plt.gca().set_facecolor('black')

    plt.tight_layout(pad=0)
    plt.show()

    return 0


def truncate_colormap(cmap, minval=0.2, maxval=1.0, n=256):
    new_cmap = LinearSegmentedColormap.from_list(
        f"trunc({cmap.name},{minval:.2f},{maxval:.2f})",
        cmap(np.linspace(minval, maxval, n))
    )
    return new_cmap


def line_plot(frame, x, y, z=None,
              cmap: Union[str, matplotlib.colors.Colormap] ='inferno',
              lw: float = 0.5):

    if z is None:
        z = np.array(x)
    # Create segments from the x and y arrays
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    # Normalize and apply colormap
    norm = Normalize(vmin=z.min(), vmax=z.max())
    lc = LineCollection(segments, cmap=cmap, norm=norm)
    lc.set_array(z)
    lc.set_linewidth(lw)
    # Set axes background to black
    frame.add_collection(lc)
    frame.autoscale()


def plot_s1d():
    # load table
    table = Table.read(S1D, hdu=1)
    # Truncate inferno to skip the darkest 20%
    trunc_cmap = truncate_colormap(INFERNO, 0.2, 1.0)

    plt.close()
    plt.figure(figsize=(32, 12), facecolor='black')  # Set figure background

    frame0 = plt.subplot2grid((3, 3), (0, 0), colspan=3)
    frame1 = plt.subplot2grid((3, 3), (1, 0), colspan=3)
    frame2 = plt.subplot2grid((3, 3), (2, 0), colspan=1)
    frame3 = plt.subplot2grid((3, 3), (2, 1), colspan=1)
    frame4 = plt.subplot2grid((3, 3), (2, 2), colspan=1)


    line_plot(frame0, table['Wave'], table['FluxAB'], cmap=trunc_cmap,
              lw=0.5)
    frame0.set(xlim=[950, 2500])

    frame1.plot(table['Wave'], table['FluxAB'], color='0.5', alpha=0.5,
                lw=0.5)
    line_plot(frame1, table['Wave'], table['FluxABTelluCorrected'],
              cmap=trunc_cmap, lw=0.5)
    frame1.set(xlim=[950, 2500])

    zoom_frames = [frame2, frame3, frame4]

    for b_it, band in enumerate(BANDS):

        wmask = (table['Wave'] > band[0]) & (table['Wave'] < band[1])
        ymin = np.nanmin(table['FluxAB'][wmask])
        ymax = np.nanmax(table['FluxAB'][wmask])
        ydiff = ymax - ymin

        zoom_frames[b_it].plot(table['Wave'], table['FluxAB'],
                               color='0.5', alpha=0.5, lw=1.0)
        line_plot(zoom_frames[b_it], table['Wave'],
                  table['FluxABTelluCorrected'],
                  cmap=trunc_cmap, lw=1.0)
        zoom_frames[b_it].set(xlim=band,
                              ylim=[ymin - 0.5* ydiff, ymax + 0.5 * ydiff])
        zoom_frames[b_it].xaxis.set_major_locator(MaxNLocator(nbins=3))

    # Turn off axes
    for frame in [frame0, frame1, frame2, frame3, frame4]:
        # Only x-ticks visible, in white
        frame.tick_params(axis='x', colors='white', labelsize=FONTSIZE)
        frame.tick_params(axis='y', which='both', left=False,
                          right=False, labelleft=False)
        frame.set_facecolor('black')
        # Axis line (spine)
        frame.spines['bottom'].set_color('white')  # x-axis line in white
        frame.spines['top'].set_color('none')
        frame.spines['left'].set_color('none')
        frame.spines['right'].set_color('none')

    plt.subplots_adjust(hspace=0.1, wspace=0.2, left=0.05, right=0.95,
                        top=0.95, bottom=0.05)

    frame0.set_title('SPIRou: Gl699', color='white', fontsize=FONTSIZE)
    plt.show()

    return 0



def plot_rdb():

    table = Table.read(RDB)

    t = Time(table['rjd'], format='mjd')
    x = t.plot_date
    y0 = table['vrad']
    yerr0 = table['svrad']
    z = t.mjd

    y1 = table['DTEMP3000']
    yerr1 = table['sDTEMP3000']

    # get RV uncertainty
    rv_val = np.nanmedian(y0)
    rv_uncertainty = np.nanmedian(yerr0)
    rv_delta_time = (np.nanmax(z) - np.nanmin(z)) / 365.254

    dtemp_val = np.nanmedian(y1)
    dtemp_uncertainty = np.nanmedian(yerr1)

    # Create colormap
    norm = plt.Normalize(vmin=np.min(z), vmax=np.max(z))

    # Truncate inferno to skip the darkest 20%
    trunc_cmap = truncate_colormap(INFERNO, 0.2, 1.0)
    colors = trunc_cmap(norm(z))

    plt.close()
    fig, frames = plt.subplots(nrows=2, figsize=(16, 16), facecolor='black',
                               sharex='all')

    # Plot error bars manually
    for xi, yi, yerri, ci in zip(x, y0, yerr0, colors):
        frames[0].errorbar(xi, yi, yerr=yerri, fmt='.',
                          color=ci, ecolor=ci, elinewidth=1, capsize=2,
                          alpha=0.5)

    for xi, yi, yerr, ci in zip(x, y1, yerr1, colors):
        frames[1].errorbar(xi, yi, yerr=yerri, fmt='.',
                          color=ci, ecolor=ci, elinewidth=1, capsize=2,
                          alpha=0.5)

    for frame in frames:
        # set formatter for date
        frame.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
        # set background black
        frame.set_facecolor('black')
        # axis labels in white
        frame.tick_params(axis='both', colors='white', labelsize=FONTSIZE)
        # remove offste notation
        frame.yaxis.get_offset_text().set_visible(False)  # hide the offset text
        frame.yaxis.set_major_formatter(mticker.ScalarFormatter(useOffset=False,
                                                                useMathText=False))
        # Axis line (spine)
        frame.spines['bottom'].set_color('white')  # x-axis line in white
        frame.spines['top'].set_color('white')
        frame.spines['left'].set_color('white')
        frame.spines['right'].set_color('white')
        # push x-axis into dates
        fig.autofmt_xdate()

        # set labels
        frame.set_xlabel('Date', color='white', fontsize=FONTSIZE)

    # set limits
    frames[0].set(ylim=[-110245, -110195])
    frames[1].set(ylim=[-5, 5])

    frames[0].set_ylabel('LBL RV [m/s]', color='white', fontsize=FONTSIZE)
    frames[1].set_ylabel('LBL DTemp [K]', color='white', fontsize=FONTSIZE)

    title0 = 'RV = {0:.2f} $\pm$ {1:.2f} m/s ({2:.1f} years)'
    targs0 = [rv_val, rv_uncertainty, rv_delta_time]
    frames[0].set_title(title0.format(*targs0), color='white', fontsize=FONTSIZE)

    title1 = 'DTemp = {0:.2f} $\pm$ {1:.2f} K'
    targs1 = [dtemp_val,  dtemp_uncertainty]
    frames[1].set_title(title1.format(*targs1), color='white', fontsize=FONTSIZE)

    plt.suptitle('SPIRou: GL699', color='white', fontsize=FONTSIZE)


    plt.savefig(os.path.join(PLOT_DIR, 'rdb_poster.png'))

    plt.show()




# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # raw image
    if 'raw' in PLOTS:
        raw_image()
    # plot s1d
    if 's1d' in PLOTS:
        plot_s1d()
    # plot rdb
    if 'rdb' in PLOTS:
        plot_rdb()


# =============================================================================
# End of code
# =============================================================================
