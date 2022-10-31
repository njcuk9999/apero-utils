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
NAME1 = 'md2.07248_py3.9.13'
NAME2 = 'md2.07248_py3.9.7'
NAME3 = 'md2.07254_py3.9.13'
NAME4 = 'md2.07254_py3.9.7'
NAME5 = 'neil@home'
# Define which reduction is the reference reduction
REF_NAME = str(NAME5)
# -----------------------------------------------------------------------------
# just add another entry here
#  i.e. paths[NAME3] = path/to/reduced/dir
paths = dict()
paths[NAME1] = '/scratch2/spirou/drs-data/minidata2_07248_py_3_9_13/red'
paths[NAME2] = '/scratch2/spirou/drs-data/minidata2_07248_py_3_9_7/red'
paths[NAME3] = '/scratch2/spirou/drs-data/minidata2_07XXX_py_3_9_13/red'
paths[NAME4] = '/scratch2/spirou/drs-data/minidata2_07XXX_py_3_9_7/red'
paths[NAME5] = '/scratch2/spirou/drs-data/minidata2_neilhome/red'
# -----------------------------------------------------------------------------
# add a color for each reduction (i.e. b, g, r, k, m, c, orange, purple)
COLORS = dict()
COLORS[NAME1] = 'b'
COLORS[NAME2] = 'r'
COLORS[NAME3] = 'g'
COLORS[NAME4] = 'purple'
COLORS[NAME5] = 'yellow'
# add a marker for each reduction (i.e. o, x, +, v, ^, d, s, .)
MARKERS = dict()
MARKERS[NAME1] = 'o'
MARKERS[NAME2] = 's'
MARKERS[NAME3] = '+'
MARKERS[NAME4] = 'x'
MARKERS[NAME5] = '^'
# markers needing facecolor
has_face = ['o', 's', '^']
# -----------------------------------------------------------------------------
# objects to consider
OBJECTS = ['GL699']
# define rv key
RV_KEY = 'RV_OBJ'


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

        if MARKERS[name] in has_face:
            pkwargs = dict(marker=MARKERS[name], markerfacecolor='None',
                           markeredgecolor=COLORS[name], ls='None')
        else:
            pkwargs = dict(marker=MARKERS[name], color=COLORS[name],
                           ls='None')

        frames[0].plot(bjds[name], rvs[name], label=name, **pkwargs)

        if name == REF_NAME:
            continue

        diff = rvs[REF_NAME] - rvs[name]
        diffrms = np.nanstd(diff)

        frames[1].plot(bjds[REF_NAME], diff, label=f'{REF_NAME}-{name}',
                       **pkwargs)

        print(f'rms of difference {REF_NAME}-{name} = {diffrms} m/s')

    frames[0].legend(loc=0)
    frames[0].set(xlabel='BJD', ylabel='RV m/s')
    plt.suptitle(f'Difference {RV_KEY} {OBJECTS}')
    frames[1].legend(loc=0)
    frames[1].set(xlabel='BJD', ylabel='$\Delta$RV m/s')
    plt.show()
    plt.close()

# =============================================================================
# End of code
# =============================================================================