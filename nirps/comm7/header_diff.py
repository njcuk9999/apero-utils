#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-11-20 at 19:19

@author: cook
"""
import os
import warnings

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import LinearStretch
from astropy.visualization import imshow_norm, ZScaleInterval
from mpl_toolkits.axes_grid1 import make_axes_locatable

# =============================================================================
# Define variables
# =============================================================================
fibers = ['A', 'B']
modes = ['HA', 'HE']
# --------------------------------------------------------------------------------------
labels = ['center', 'left', 'right', 'red', 'blue']
boxsize = 5
function = np.nanmedian
func_name = 'Med'
# --------------------------------------------------------------------------------------
PLOT = True
PLOT_PATH = '/nirps_raw/nirps/misc/comm7'
# --------------------------------------------------------------------------------------
RAW_DIR_HE = '/nirps_raw/nirps/raw-data/nirps_he/'
RAW_DIR_HA = '/nirps_raw/nirps/raw-data/nirps_ha/'
# --------------------------------------------------------------------------------------
# switch files key = (mode, fiber)
CASE = 1
old_files = dict()
new_files = dict()
xpos_all = dict()
ypos_all = dict()

if CASE == 1:
    # NIRPS HE Fiber A
    old_files[('HE', 'A')] = RAW_DIR_HE + '2022-11-17/NIRPS_2022-11-17T22_15_02_461.fits'
    new_files[('HE', 'A')] = RAW_DIR_HE + '2022-11-19/NIRPS_2022-11-19T23_32_10_890.fits'
    # NIRPS HE Fiber B
    old_files[('HE', 'B')] = RAW_DIR_HE + '2022-11-17/NIRPS_2022-11-17T22_16_42_792.fits'
    new_files[('HE', 'B')] = RAW_DIR_HE + '2022-11-19/NIRPS_2022-11-19T23_33_51_221.fits'
    # NIRPS HA Fiber A
    old_files[('HA', 'A')] = RAW_DIR_HA + '2022-11-17/NIRPS_2022-11-17T15_22_57_043.fits'
    new_files[('HA', 'A')] = RAW_DIR_HA + '2022-11-19/NIRPS_2022-11-19T23_28_44_649.fits'
    # NIRPS HA Fiber B
    old_files[('HA', 'B')] = RAW_DIR_HA + '2022-11-17/NIRPS_2022-11-17T15_24_37_374.fits'
    new_files[('HA', 'B')] = RAW_DIR_HA + '2022-11-19/NIRPS_2022-11-19T23_30_24_979.fits'
elif CASE == 2:
    # NIRPS HE Fiber A
    old_files[('HE', 'A')] = RAW_DIR_HE + '2022-11-19/NIRPS_2022-11-19T23_32_10_890.fits'
    new_files[('HE', 'A')] = RAW_DIR_HE + '2022-11-21/NIRPS_2022-11-21T16_46_59_021.fits'
    # NIRPS HE Fiber B
    old_files[('HE', 'B')] = RAW_DIR_HE + '2022-11-19/NIRPS_2022-11-19T23_33_51_221.fits'
    new_files[('HE', 'B')] = RAW_DIR_HE + '2022-11-21/NIRPS_2022-11-21T16_48_39_353.fits'
    # NIRPS HA Fiber A
    old_files[('HA', 'A')] = RAW_DIR_HA + '2022-11-19/NIRPS_2022-11-19T23_28_44_649.fits'
    new_files[('HA', 'A')] = RAW_DIR_HA + '2022-11-21/NIRPS_2022-11-21T16_43_38_361.fits'
    # NIRPS HA Fiber B
    old_files[('HA', 'B')] = RAW_DIR_HA + '2022-11-19/NIRPS_2022-11-19T23_30_24_979.fits'
    new_files[('HA', 'B')] = RAW_DIR_HA + '2022-11-21/NIRPS_2022-11-21T16_45_18_691.fits'
else:
    raise ValueError(f'Case = {CASE} not supported')


# =============================================================================
# Define functions
# =============================================================================

# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    storage = dict()
    comments = dict()

    for mode in modes:
        for fiber in fibers:
            old_file = old_files[(mode, fiber)]
            new_file = new_files[(mode, fiber)]


            # load file headers
            old_hdr = fits.getheader(old_file)
            new_hdr = fits.getheader(new_file)

            # storage
            different = dict()
            comments[(mode, fiber)] = dict()

            # check keys
            for key in new_hdr:
                if new_hdr[key] != old_hdr[key]:
                    try:
                        ratio = new_hdr[key] / old_hdr[key]

                        if 0.9 < ratio < 1.1:
                            continue

                        different[key] = (new_hdr[key], old_hdr[key], ratio)
                        comments[(mode, fiber)][key] = new_hdr.comments[key]
                    except:
                        continue

            storage[(mode, fiber)] = different


    for name in storage:
        mode, fiber = name
        print('*' * 50)
        print(f'NIRPS {mode} Fiber {fiber}')
        print('*' * 50)

        for key in storage[name]:
            new, old, ratio = storage[name][key]

            print(f'{key} old={old} new={new} ratio={ratio}\t\t\t// {comments[name][key]}')






# =============================================================================
# End of code
# =============================================================================
