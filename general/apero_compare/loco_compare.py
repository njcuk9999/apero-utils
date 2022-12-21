#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-10-28 at 11:42

@author: cook
"""
import os

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits


# =============================================================================
# Define variables
# =============================================================================
import os
from typing import List
import itertools
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits


# =============================================================================
# Define variables
# =============================================================================
NAME1 = 'spirou@rali'
NAME2 = 'cook@jupiter'
NAME3 = 'cook@nb19'
# Define which reduction is the reference reduction
REF_NAME = str(NAME1)
# -----------------------------------------------------------------------------
# just add another entry here
#  i.e. paths[NAME3] = path/to/reduced/dir
paths = dict()
paths[NAME1] = '/scratch2/spirou/drs-data/minidata2_2022-12-15_rali/red'
paths[NAME2] = '/scratch2/spirou/drs-data/minidata2_07XXX/reduced'
paths[NAME3] = '/scratch2/spirou/drs-data/minidata2_2022-12-15_nb19/red'
# -----------------------------------------------------------------------------
# add a color for each reduction (i.e. b, g, r, k, m, c, orange, purple)
COLORS = dict()
COLORS[NAME1] = 'b'
COLORS[NAME2] = 'r'
COLORS[NAME3] = 'g'
# COLORS[NAME4] = 'orange'
# add a marker for each reduction (i.e. o, x, +, v, ^, d, s, .)
MARKERS = dict()
MARKERS[NAME1] = 'o'
MARKERS[NAME2] = 's'
MARKERS[NAME3] = '+'
# MARKERS[NAME4] = '^'
# markers needing facecolor
has_face = ['o', 's', '^', 'd', 'v']
# -----------------------------------------------------------------------------
# objects to consider
OBJECTS = ['GL699']
# file suffix to identify file
FILE_SUFFIX = 'loco_AB.fits'
# plotting parameters
N_X_POINTS = 4
FIGSHOW = True
FIGNAME = None

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




def loco_plot(lfiles):
    # get the number of plots if we want a square grid
    indices = int(np.ceil(np.sqrt(N_X_POINTS)))
    # get the positions of each plot
    graph_pos = list(itertools.product(range(indices), range(indices)))
    # get the number of columns and rows in plot
    cols, rows = indices, int(np.ceil(N_X_POINTS / indices))
    # deal with fibers
    if lfiles[REF_NAME].endswith('AB.fits'):
        fibers = dict(A=0, B=1)
    else:
        fibers = dict(C=None)
    # get the reference data
    data_ref = fits.getdata(lfiles[REF_NAME])
    # get the x pixel positions to plot at
    xpixels = np.linspace(0, data_ref.shape[1] - 1, N_X_POINTS, dtype=int)

    # -------------------------------------------------------------------------
    # loop around each fiber
    for fiber in fibers:
        # get the frame to plot on
        plt.close()
        fig, frames = plt.subplots(ncols=cols, nrows=rows,
                                   figsize=(5 * cols, 5 * rows))
        # loop round each file
        for name_cmp in lfiles:
            # don't compare ref_name to ref_name
            if name_cmp == REF_NAME:
                continue
            print(f'Plotting {REF_NAME} - {name_cmp}')
            # get the comparison data
            data_cmp = fits.getdata(lfiles[name_cmp])
            if fibers[fiber] == 0:
                pdata1 = data_ref[0::2]
                pdata2 = data_cmp[0::2]
            elif fibers[fiber] == 1:
                pdata1 = data_ref[1::2]
                pdata2 = data_cmp[1::2]
            else:
                pdata1 = data_ref
                pdata2 = data_cmp

            if MARKERS[name_cmp] in has_face:
                pkwargs = dict(marker=MARKERS[name_cmp], markerfacecolor='None',
                               markeredgecolor=COLORS[name_cmp], ls='None')
            else:
                pkwargs = dict(marker=MARKERS[name_cmp], color=COLORS[name_cmp],
                               ls='None')

            # loop around x pixel positions
            for it, xpixel in enumerate(xpixels):
                # deal with fibers (AB means we have to get A and B)
                # C is easier
                frame = frames[graph_pos[it]]
                diff = pdata1[:, xpixel]-pdata2[:, xpixel]
                rmsdiff = np.nanstd(diff)
                # plot fiber
                frame.plot(np.arange(pdata1.shape[0]), diff,
                           label=f'{REF_NAME}-{name_cmp}', **pkwargs)

                print(f'RMS[{fiber}][{xpixel}]: {rmsdiff}')

        # loop around x pixel positions
        for it, xpixel in enumerate(xpixels):
            # get the frame to plot on
            frame = frames[graph_pos[it]]
            frame.set(xlabel='Order', ylabel=r'$\Delta$ Y pixel',
                      title=f'X pixel = {xpixel}')
            frame.legend(loc=0)

        for it in range(N_X_POINTS, frames.size):
            frame = frames[graph_pos[it]]
            frame.axis('off')

        plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05,
                            top=0.95, hspace=0.2, wspace=0.2)
        plt.suptitle(lfiles[REF_NAME] + f' FIBER = {fiber}')
        if FIGNAME is not None:
            print(f'Saved to: {FIGNAME}')
            plt.savefig(FIGNAME)
        if FIGSHOW:
            print(f'Showing plot {FIGNAME}')
            plt.show()
        plt.close()


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

        if len(files[name]) == 0:
            print(f'No files in {path}')
    # -------------------------------------------------------------------------
    # only deal with files with e2dsff_tcorr_AB.fits
    # we take the first "name" as the reference
    basenames = []
    for filename in files[REF_NAME]:
        if filename.endswith(FILE_SUFFIX):
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
    # plot loco plot
    # loop around valid files
    for valid_basename in valid_basenames:

        valid_files = dict()
        for name in paths:
            valid_files[name] = os.path.join(paths[name], valid_basename)

        # plot loco plot
        loco_plot(valid_files)




# =============================================================================
# End of code
# =============================================================================
