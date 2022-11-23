#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-11-20 at 19:19

@author: cook
"""
import os
import glob

import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from astropy.io import fits
from tqdm import tqdm


# =============================================================================
# Define variables
# =============================================================================
# raw directories
RAW_DIR_HE = '/nirps_raw/nirps/raw-data/nirps_he/'
RAW_DIR_HA = '/nirps_raw/nirps/raw-data/nirps_ha/'

# observation directory
OBS_DIR = '2022-11-21'

# define start and end time of test (only use files in this range)
TIME_START = Time('2022-11-21 23:40:00')
TIME_END = Time('2022-11-22 00:45:00')
# define fibers and modes
fibers = ['A', 'B']
modes = ['HA', 'HE']
# header keys
KW_RAW_DPRTYPE = 'HIERARCH ESO DPR TYPE'
KW_ND_FILT1 = 'HIERARCH ESO INS FILT1 NAME'
KW_ND_FILT2 = 'HIERARCH ESO INS FILT2 NAME'

# =============================================================================
# Define functions
# =============================================================================
def convert_nd_to_value(filename, input_value):
    try:
        return float(input_value[2:])
    except:
        ValueError(f'ND cannot be converted ND={input_value} filename={filename}')


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # storage of data
    nd_arr = []
    mode_arr = []
    fiber_arr = []
    per90_arr = []
    # ----------------------------------------------------------------------
    # loop around modes and all files in each raw obs dir
    for mode in modes:
        print('=' * 50 + f'NIRPS_{mode}'  + '=' * 50)

        # get all files for this obs_dir
        if mode == 'HA':
            files = glob.glob(os.path.join(RAW_DIR_HA, OBS_DIR, '*.fits'))
        else:
            files = glob.glob(os.path.join(RAW_DIR_HE, OBS_DIR, '*.fits'))
        # ------------------------------------------------------------------
        # get data and header for each fiber
        for filename in tqdm(files):
            # load header
            header = fits.getheader(filename)
            # ------------------------------------------------------------------
            raw_dprtype = header[KW_RAW_DPRTYPE]
            # only interested in order def files
            if raw_dprtype not in ['ORDERDEF,LAMP,DARK', 'ORDERDEF,DARK,LAMP']:
                continue
            # ------------------------------------------------------------------
            # get date
            date = Time(header['DATE'], format='fits')
            # skip files not in test range
            if date < TIME_START or date > TIME_END:
                continue
            # ------------------------------------------------------------------
            # we have fiber A
            if raw_dprtype == 'ORDERDEF,LAMP,DARK':
                data = fits.getdata(filename)
                fiber_arr.append('A')
                nd_val = convert_nd_to_value(filename, header[KW_ND_FILT1])
                nd_arr.append(nd_val)
                mode_arr.append(mode)
                per90_arr.append(np.nanpercentile(data, 90))
            # we have fiber B
            elif raw_dprtype == 'ORDERDEF,DARK,LAMP':
                data = fits.getdata(filename)
                fiber_arr.append('B')
                nd_val = convert_nd_to_value(filename, header[KW_ND_FILT2])
                nd_arr.append(nd_val)
                mode_arr.append(mode)
                per90_arr.append(np.nanpercentile(data, 90))
    # ----------------------------------------------------------------------
    # convert to numpy arrays
    nd_arr = np.array(nd_arr)
    mode_arr = np.array(mode_arr)
    fiber_arr = np.array(fiber_arr)
    per90_arr = np.array(per90_arr)
    # ----------------------------------------------------------------------
    # now we can plot
    fig, frame = plt.subplots(ncols=1, nrows=1)

    markers = ['x', '+', 'o', '*']
    edgecolor = ['r', 'g', 'b', 'purple']
    facecolor = ['r', 'g', 'None', 'None']
    it = 0
    for mode in modes:
        for fiber in fibers:
            # set up the mask for this combination
            mask = (mode_arr == mode) & (fiber_arr == fiber)

            x_arr, y_arr = nd_arr[mask], per90_arr[mask]

            sortmask = np.argsort(x_arr)
            x_arr, y_arr = x_arr[sortmask], y_arr[sortmask]
            # plot this mask
            frame.plot(x_arr, y_arr, label=f'{mode}-{fiber}',
                       markeredgecolor=edgecolor[it], markerfacecolor=facecolor[it],
                       marker=markers[it], color=edgecolor[it])
            it += 1

    frame.set(xlabel='ND value', ylabel='90% percentile flux (full image)')
    frame.legend(loc=0)
    plt.show()
    plt.close()

# =============================================================================
# End of code
# =============================================================================
