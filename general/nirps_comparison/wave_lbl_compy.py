#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-06-10 at 12:31

@author: cook
"""
import os

import numpy as np
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from astropy.io import fits
from astropy.table import Table
from astropy import constants
from astropy import units as uu
from tqdm import tqdm
import warnings
from scipy.spatial import cKDTree


# =============================================================================
# Define variables
# =============================================================================
# Set global font size
mpl.rcParams.update({'font.size': 24})

BASENAME = 'NIRPS.2025-06-07T07:04:20.565'

ORDER_NUM = 60

speed_of_light_ms = constants.c.value

# Set global font size
mpl.rcParams.update({'font.size': 24})
# Set background color in plots
PLOT_BACKGROUND_COLOR = '#FEFDE1'
PLOT_BACKGROUND_COLOR2 = '#E1F0FE'

# set plot path
PLOT_PATH = '/data/cook/nirps_comp/plots'

MJD_KEY = 'MJD-OBS'

ORDER_DEF_CACHE = None


# =============================================================================
# Define classes
# =============================================================================
class File:
    def __init__(self, name):
        self.path = ''
        self.obsdir = None
        self.basename = ''
        self.prefix = ''
        self.suffix = ''
        self.format = None

    def get_abspath(self):
        basename = self.prefix + self.basename + self.suffix
        if self.obsdir is None:
            abspath = os.path.join(self.path, basename)
        else:
            abspath = os.path.join(self.path, self.obsdir, basename)
        if not os.path.exists(abspath):
            raise ValueError(f'File {abspath} does not exist.')
        return abspath

    def load_fits(self):
        abspath = self.get_abspath()

        if not os.path.exists(abspath):
            raise FileNotFoundError(abspath)

        return Table.read(abspath)

# =============================================================================
# Define functions
# =============================================================================
def match_lines(lines1, lines2, wid_lines1, wid_lines2):
    """
    Match lines between two sets by nearest wavelength within a max difference.
    Returns list of matched pairs (λ1, λ2)
    """
    tree = cKDTree(np.array(lines2).reshape(-1, 1))
    matches = []

    for idx1, lambda1 in enumerate(lines1):
        # get the distance all lines1 are from lines2
        dist, all_idx2s = tree.query([[lambda1]], k=1)
        # first one is the closest
        idx2 = all_idx2s[0]
        # take the max size of this line and
        max_size = np.max([wid_lines1[idx1], wid_lines2[idx2]])
        # if the line is inside the limits of line1 and line2 then keep it
        if dist[0] <= max_size / 2:
            lambda2 = lines2[idx2]
            matches.append((lambda1, lambda2))

    return np.array(matches)


def apero2eso_ord(order_num: int, fmt='apero') -> int:
    """Convert APERO order number to ESO order number."""

    # load csv
    order_def = Table.read('order_def_apero_eso.csv')

    if fmt == 'apero':
        order_dict = dict(zip(order_def['APERO_PYTHON'],
                              order_def['ESO_PYTHON']))
    elif fmt == 'eso':
        order_dict = dict(zip(order_def['ESO_PYTHON'],
                              order_def['APERO_PYTHON']))
    else:
        raise ValueError(f'Format {fmt} not recognized. Use "apero" or "eso".')

    if order_num not in order_dict:
        raise ValueError(f'Order number {order_num} not found in order definition.')

    if np.isnan(order_dict[order_num]):
        raise ValueError(f'Order number {order_num} is not defined in order definition.')

    return int(order_dict[order_num])


def get_echelle_ord(order_num: int, fmt='apero') -> int:

    global ORDER_DEF_CACHE
    # load csv
    if ORDER_DEF_CACHE is not None:
        order_def = ORDER_DEF_CACHE
    else:
        order_def = Table.read('order_def_apero_eso.csv')
        ORDER_DEF_CACHE = order_def

    if fmt == 'apero':
        order_dict = dict(zip(order_def['APERO_PYTHON'],
                              order_def['ECHELLE_ORDER']))
    elif fmt == 'eso':
        order_dict = dict(zip(order_def['ESO_PYTHON'],
                              order_def['ECHELLE_ORDER']))
    else:
        raise ValueError(f'Format {fmt} not recognized. Use "apero" or "eso".')

    if order_num not in order_dict:
        raise ValueError(f'Order number {order_num} not found in order definition.')

    if np.isnan(order_dict[order_num]):
        raise ValueError(f'Order number {order_num} is not defined in order definition.')

    return int(order_dict[order_num])


# =============================================================================
# Define files
# =============================================================================
apero = File('apero')
apero.path = '/scratch2/nirps/misc/nirps_comp/lbl/apero/lblrv/GL699_GL699'
apero.basename = BASENAME
apero.suffix = 't_GL699_GL699_lbl.fits'
apero.format = 'table'

eso = File('eso')
eso.path = '/scratch2/nirps/misc/nirps_comp/lbl/eso/lblrv/GL699_GL699'
eso.basename = BASENAME
eso.prefix = 'r.'
eso.suffix = '_S2D_TELL_CORR_A_GL699_GL699_lbl.fits'
eso.format = 'table'

# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # load tables
    apero_table = apero.load_fits()
    eso_table = eso.load_fits()

    apero_chunk_cent = np.nanmean([apero_table['WAVE_END'], apero_table['WAVE_START']], axis=0)
    apero_chunk_wid = apero_table['WAVE_END'] - apero_table['WAVE_START']
    eso_chunk_cent = np.nanmean([eso_table['WAVE_END'], eso_table['WAVE_START']], axis=0)
    eso_chunk_wid = eso_table['WAVE_END'] - eso_table['WAVE_START']

    # match lines between catalgoues
    matches = match_lines(apero_chunk_cent, eso_chunk_cent,
                          apero_chunk_wid / 10, eso_chunk_wid / 10)
    # get the difference in wavelength
    waves_diff = matches[:, 1] - matches[:, 0]
    # get the mean wavelength between the two
    waves_mean = (matches[:, 1] + matches[:, 0]) / 2

    # get the mean wave diff (bulk shift)
    mean_wave, std_wave = np.nanmean(waves_diff), np.nanstd(waves_diff)
    # get the difference in m/s
    rv_diff = speed_of_light_ms * waves_diff / waves_mean
    # get the mean rv diff (bulk shift)
    mean_rv, std_rv = np.nanmean(rv_diff), np.nanstd(rv_diff)

    # -------------------------------------------------------------------------
    # plot
    fig, frames = plt.subplots(nrows=1, ncols=2)

    # plot the histogram of the pixel difference
    frames[0].hist(waves_diff, bins=1000, color='orange')
    frames[0].set(xlabel='Wavelength center difference [nm]',
                  title=f'Bulk shift = {mean_wave:.3e}$\pm${std_wave:.3e} nm',
                  xlim=[-5*std_wave, 5*std_wave])

    frames[1].hist(rv_diff, bins=1000, color='blue')
    frames[1].set(xlabel='RV difference [m/s]',
                  title=f'Bulk shift = {mean_rv:.2f}$\pm${std_rv:.2f} m/s',
                  xlim=[-5*std_rv, 5*std_rv])


    plt.savefig(os.path.join(PLOT_PATH, 'wave_lbl_comp.png'))
    plt.show()

    # find orders for each match
    order_matches = np.full(matches.shape, np.nan)

    # loop around rows
    for row in tqdm(range(len(order_matches))):

        # find position in apero cents
        apero_pos = np.where(apero_chunk_cent == matches[row, 0])[0][0]
        eso_pos = np.where(eso_chunk_cent == matches[row, 1])[0][0]

        # get order positions for both
        order_matches[row, 0] = get_echelle_ord(apero_table['ORDER'][apero_pos], 'apero')
        order_matches[row, 1] = get_echelle_ord(eso_table['ORDER'][eso_pos], 'eso')

    # -------------------------------------------------------------------------
    # make a histogram of every order
    min_order = np.nanmin(order_matches)
    max_order = np.nanmax(order_matches)

    orders = np.arange(min_order, max_order + 1).astype(int)

    rv_diff_image = np.zeros([len(orders), 1001])


    for it, order_num in enumerate(orders):

        # mask orders
        order_mask = order_matches == order_num

        rv_mask = np.mean(order_mask, axis=1) == 1.0

        hist, _ = np.histogram(rv_diff[rv_mask], range=[-1000, 1000])

        if len(hist) == 1001:
            rv_diff_image[it] = hist




# =============================================================================
# End of code
# =============================================================================
