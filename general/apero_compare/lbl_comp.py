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
from astropy.table import Table

# =============================================================================
# Define variables
# =============================================================================
NAME1 = 'spirou@rali'
NAME2 = 'cook@jupiter'
NAME3 = 'cook@nb19'
NAME4 = 'cook@home'
# Define which reduction is the reference reduction
REF_NAME = str(NAME1)
# -----------------------------------------------------------------------------
# just add another entry here
#  i.e. paths[NAME3] = path/to/reduced/dir
paths = dict()
paths[NAME1] = '/scratch3/lbl/data/minidata_comp/minidata2_2022-12-15_rali/lblrdb/'
paths[NAME2] = '/scratch3/lbl/data/minidata_comp/minidata2_07XXX/lblrdb/'
paths[NAME3] = '/scratch3/lbl/data/minidata_comp/minidata2_2022-12-15_nb19/lblrdb/'
# paths[NAME4] = '/scratch3/lbl/data/minidata_comp/minidata2_2022-12-15_home/lblrdb/'
# -----------------------------------------------------------------------------
# add a color for each reduction (i.e. b, g, r, k, m, c, orange, purple)
COLORS = dict()
COLORS[NAME1] = 'b'
COLORS[NAME2] = 'r'
COLORS[NAME3] = 'g'
COLORS[NAME4] = 'orange'
# add a marker for each reduction (i.e. o, x, +, v, ^, d, s, .)
MARKERS = dict()
MARKERS[NAME1] = 'o'
MARKERS[NAME2] = 's'
MARKERS[NAME3] = '+'
MARKERS[NAME4] = '^'
# markers needing facecolor
has_face = ['o', 's', '^', 'd', 'v']
# -----------------------------------------------------------------------------
# objects to consider
OBJECTS = ['GL699']
# define rv key
TIME_KEY = 'rjd'
RV_KEY = 'vrad'
ERV_KEY = 'svrad'


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
    objnames = []
    for filename in files[REF_NAME]:
        if filename.endswith('.fits'):
            basename = filename.split(paths[REF_NAME])[-1]

            while basename.startswith(os.sep):
                basename = basename[len(os.sep):]

            basenames.append(basename)
            objname = basename.replace('.fits', '').split('_')[-2]
            objnames.append(objname)
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
    times = dict()
    rvs = dict()
    ervs = dict()
    obsdirs = dict()

    # loop around reductions
    for it, name in enumerate(paths):
        # define a list per reduction to store bjds and rvs
        times[name] = []
        rvs[name] = []
        ervs[name] = []
        obsdirs[name] = []
        # loop around valid basenames
        for jt, valid_basename in enumerate(valid_basenames):
            try:
                file_n = os.path.join(paths[name], valid_basename)
                table = Table.read(file_n, hdu=9)
                hdr = fits.getheader(file_n, hdu=0)
            except Exception as e:
                print(f'Cannot load {valid_basename} for {name} -> skipping')
                print(f'\tError {type(e)}: {str(e)}')
                times[name].append(np.nan)
                ervs[name].append(np.nan)
                rvs[name].append(np.nan)
                continue
            # filter out unwanted object types
            if objnames[jt] not in OBJECTS:
                times[name].append(np.nan)
                ervs[name].append(np.nan)
                rvs[name].append(np.nan)
            else:
                times[name] = table[TIME_KEY]
                rvs[name] = table[RV_KEY]
                ervs[name] = table[ERV_KEY]
    # -------------------------------------------------------------------------
    # keep track of nans
    nanmask = np.zeros(len(times[REF_NAME]), dtype=bool)
    # convert to numpy arrays and count nanrows
    for name in paths:
        times[name] = np.array(times[name])
        rvs[name] = np.array(rvs[name])
        ervs[name] = np.array(ervs[name])
        # we need a combined mask between reductions
        nanmask |= ~np.isfinite(times[name])
        nanmask |= ~np.isfinite(rvs[name])
        nanmask |= ~np.isfinite(ervs[name])

    # -------------------------------------------------------------------------
    # plot
    fig, frames = plt.subplots(ncols=1, nrows=2)

    for name in paths:

        xvector = times


        if MARKERS[name] in has_face:
            pkwargs = dict(marker=MARKERS[name], markerfacecolor='None',
                           markeredgecolor=COLORS[name], ls='None')
        else:
            pkwargs = dict(marker=MARKERS[name], color=COLORS[name],
                           ls='None')

        frames[0].errorbar(xvector[name], rvs[name], yerr=ervs[name],
                           label=name, **pkwargs)

        if name == REF_NAME:
            continue

        diff = rvs[REF_NAME] - rvs[name]
        ediff = np.sqrt(ervs[REF_NAME]**2 + ervs[name]**2)
        diffrms = np.nanstd(diff)

        frames[1].plot(xvector[REF_NAME], diff,
                       label=f'{REF_NAME}-{name}', **pkwargs)

        print(f'rms of difference {REF_NAME}-{name} = {diffrms} m/s')

    frames[0].legend(loc=0)

    frames[0].set(xlabel=TIME_KEY, ylabel='RV m/s')
    plt.suptitle(f'Difference {RV_KEY} {OBJECTS}')
    frames[1].legend(loc=0)
    frames[1].set(xlabel=TIME_KEY, ylabel='$\Delta$RV m/s')
    plt.show()
    plt.close()

# =============================================================================
# End of code
# =============================================================================
