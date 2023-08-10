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
name1 = 'spirou@rali'
name2 = 'cook@jupiter'
name3 = 'cook@nb19'
name4 = 'spirou@maestria'
name5 = 'newworlds'
name6 = 'lam'
name7 = 'cfht'


names = [name1, name2, name3, name4, name5, name6, name7]

# This is a hack but just to test without certain points
REJECT_DATE_STARTS = [59063.7786]
REJECT_DATE_ENDS = [59063.7790]

# Define which reduction is the reference reduction
REF_NAME = str(name1)
# -----------------------------------------------------------------------------
# just add another entry here
#  i.e. paths[NAME3] = path/to/reduced/dir
outpaths = dict()
outpaths[name1] = '/scratch3/lbl/data/minidata_comp/spirou_minidata2_07286_rali/lblrdb/'
outpaths[name2] = '/scratch3/lbl/data/minidata_comp/spirou_minidata2_07286_jupiter/lblrdb/'
outpaths[name3] = '/scratch3/lbl/data/minidata_comp/spirou_minidata2_07286_nb19/lblrdb/'
outpaths[name4] = '/scratch3/lbl/data/minidata_comp/spirou_minidata2_07286_maestria/lblrdb/'
outpaths[name5] = '/scratch3/lbl/data/minidata_comp/spirou_minidata2_07286_newworld/lblrdb'
outpaths[name6] = '/scratch3/lbl/data/minidata_comp/spirou_minidata2_07286_lam/lblrdb/'
outpaths[name7] = '/scratch3/lbl/data/minidata_comp/spirou_minidata2_07286_cfht/lblrdb/'

paths = outpaths

# -----------------------------------------------------------------------------
# add a color for each reduction (i.e. b, g, r, k, m, c, orange, purple)
COLORS = dict()
COLORS[name1] = 'b'
COLORS[name2] = 'r'
COLORS[name3] = 'g'
COLORS[name4] = 'orange'
COLORS[name5] = 'purple'
COLORS[name6] = 'k'
COLORS[name7] = 'm'
# add a marker for each reduction (i.e. o, x, +, v, ^, d, s, .)
MARKERS = dict()
MARKERS[name1] = 'o'
MARKERS[name2] = 's'
MARKERS[name3] = '+'
MARKERS[name4] = '^'
MARKERS[name5] = 'x'
MARKERS[name6] = 'v'
MARKERS[name7] = 'd'
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
    for name in names:
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
        for name in names:
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
    for it, name in enumerate(names):
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
    for name in names:
        times[name] = np.array(times[name])
        rvs[name] = np.array(rvs[name])
        ervs[name] = np.array(ervs[name])
        # we need a combined mask between reductions
        nanmask |= ~np.isfinite(times[name])
        nanmask |= ~np.isfinite(rvs[name])
        nanmask |= ~np.isfinite(ervs[name])
    # -------------------------------------------------------------------------
    # reject points
    if len(REJECT_DATE_ENDS) > 0:

        reject_mask = np.zeros_like(times[REF_NAME], dtype=bool)
        # convert to numpy arrays and count nanrows
        for name in names:
            for point in range(len(REJECT_DATE_ENDS)):
                point_mask = times[name] > REJECT_DATE_STARTS[point]
                point_mask &= times[name] < REJECT_DATE_ENDS[point]

                reject_mask |= point_mask
        # apply mask
        for name in names:
            times[name] = times[name][~reject_mask]
            rvs[name] = rvs[name][~reject_mask]
            ervs[name] = ervs[name][~reject_mask]
    # -------------------------------------------------------------------------
    # plot
    fig, frames = plt.subplots(ncols=1, nrows=3)

    for name in names:

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

        diff2 = ((rvs[REF_NAME] - np.nanmedian(rvs[REF_NAME])) -
                 (rvs[name] - np.nanmedian(rvs[name])))


        frames[1].plot(xvector[REF_NAME], diff,
                       label=f'{REF_NAME}-{name}', **pkwargs)

        frames[2].plot(xvector[REF_NAME], diff2,
                       label=f'({REF_NAME}-med)-({name}-med)', **pkwargs)

        print(f'rms of difference {REF_NAME}-{name} = {diffrms} m/s')

    frames[0].legend(loc=0)
    frames[0].ticklabel_format(useOffset=False)

    frames[0].set(xlabel=TIME_KEY, ylabel='RV m/s')
    plt.suptitle(f'Difference {RV_KEY} {OBJECTS}')
    frames[1].legend(loc=0)
    frames[1].set(xlabel=TIME_KEY, ylabel='$\Delta$RV m/s')
    frames[1].ticklabel_format(useOffset=False)

    frames[2].legend(loc=0)
    frames[2].set(xlabel=TIME_KEY, ylabel='$\Delta$RV m/s')
    frames[2].ticklabel_format(useOffset=False)


    plt.show()
    plt.close()

# =============================================================================
# End of code
# =============================================================================
