#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-10-27 at 17:41

@author: cook
"""
import os
from typing import List

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits


# =============================================================================
# Define variables
# =============================================================================
NAME1 = 'md2.Udem.v0.7.254'
NAME2 = 'md2.Neil.Home'
# Define which reduction is the reference reduction
REF_NAME = str(NAME2)
# -----------------------------------------------------------------------------
# just add another entry here
#  i.e. paths[NAME3] = path/to/reduced/dir
paths = dict()
paths[NAME1] = '/scratch2/spirou/drs-data/minidata2_07XXX_extmem/reduced'
paths[NAME2] = '/scratch2/spirou/drs-data/minidata2_neilhome/red'
# -----------------------------------------------------------------------------
# add a color for each reduction (i.e. b, g, r, k, m, c, orange, purple)
colors = dict()
colors[NAME1] = 'b'
colors[NAME2] = 'g'
# add a marker for each reduction (i.e. o, x, +, v, ^, d, s, .)
markers = dict()
markers[NAME1] = 'o'
markers[NAME2] = 'x'
# -----------------------------------------------------------------------------
# objects to consider
OBJECTS = ['GL699']


# =============================================================================
# Define functions
# =============================================================================
def get_files(path: str) -> List[str]:
    fits_files = []
    for root, dir, files in os.walk(path):
        for filename in files:
            if filename.endswith('.fits'):
                fits_files.append(os.path.join(root, filename))
    return fits_files


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # -------------------------------------------------------------------------


    # get all fits files
    files = dict()
    for name in paths:
        path = paths[name]
        files[name] = get_files(path)
    # -------------------------------------------------------------------------
    # only deal with files with e2dsff_tcorr_AB.fits
    # we take the first "name" as the reference
    basenames = []
    for filename in files[REF_NAME]:
        if filename.endswith('e2dsff_tcorr_AB_ccf_gl699_neg.fits_AB.fits'):
            basename = filename.split(paths[REF_NAME])[-1]

            while basename.startswith(os.sep):
                basename = basename[len(os.sep):]

            basenames.append(basename)
    # -------------------------------------------------------------------------
    # need to only consider files that appear in all reductions
    valid_basenames = []

    for basename in basenames:
        cond = True
        for name in paths:
            cond &= os.path.exists(os.path.join(paths[name], basename))
        if cond:
            valid_basenames.append(basename)

    # -------------------------------------------------------------------------
    # ge bjd vs rv for all objects
    bjds = dict()
    rvs = dict()

    # loop around reductions
    for it, name in enumerate(paths):
        # define a list per reduction to store bjds and rvs
        bjds[name] = []
        rvs[name] = []
        # loop around valid basenames
        for valid_basename in valid_basenames:
            try:
                file_n = os.path.join(paths[name], valid_basename)
                hdr = fits.getheader(file_n)
            except Exception as e:
                print(f'Cannot load {valid_basename} for {name} -> skipping')
                print(f'\tError {type(e)}: {str(e)}')
                bjds[name].append(np.nan)
                rvs[name].append(np.nan)
                continue
            # filter out unwanted object types
            if hdr['OBJECT'].upper() not in OBJECTS:
                bjds[name].append(np.nan)
                rvs[name].append(np.nan)
            else:
                bjds[name].append(hdr['BJD'])
                rvs[name].append(hdr['RV_OBJ'] * 1000)

    # -------------------------------------------------------------------------
    # keep track of nans
    nanmask = np.zeros(len(valid_basenames), dtype=bool)
    # convert to numpy arrays and count nanrows
    for name in paths:
        bjds[name] = np.array(bjds[name])
        rvs[name] = np.array(rvs[name])
        # we need a combined mask between reducions
        nanmask |= ~np.isfinite(bjds[name])
        nanmask |= ~np.isfinite(rvs[name])

    # -------------------------------------------------------------------------
    # plot
    fig, frames = plt.subplots(ncols=1, nrows=2)

    for name in paths:
        frames[0].scatter(bjds[name], rvs[name], label=name,
                          color=colors[name], marker=markers[name])

    frames[0].legend(loc=0)
    frames[0].set(xlabel='BJD', ylabel='RV m/s')

    for name in paths:
        if name == REF_NAME:
            continue
        frames[1].scatter(bjds[REF_NAME], rvs[REF_NAME] - rvs[name],
                          label=f'{REF_NAME}-{name}',
                          color=colors[name], marker=markers[name])
    frames[1].legend(loc=0)
    frames[1].set(xlabel='BJD', ylabel='$\Delta$RV m/s')
    plt.show()
    plt.close()

# =============================================================================
# End of code
# =============================================================================
