#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-11-15 at 13:57

@author: cook
"""
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

# =============================================================================
# Define variables
# =============================================================================
NAME1 = 'NW'
NAME2 = 'cook@jupiter'
NAME3 = 'cook@nb19'
NAME4 = 'spirou@rali'
# Define which reduction is the reference reduction
REF_NAME = str(NAME1)
# -----------------------------------------------------------------------------
# just add another entry here
#  i.e. paths[NAME3] = path/to/reduced/dir
paths = dict()
paths[NAME1] = '/scratch2/spirou/misc/miniruns_v07261/claire_nw.txt'
paths[NAME2] = '/scratch2/spirou/misc/miniruns_v07261/cook_jupiter.txt'
# paths[NAME3] = '/scratch2/spirou/misc/miniruns_v07261/cook_nb19.txt'
paths[NAME4] = '/scratch2/spirou/misc/miniruns_v07261/spirou_rali.txt'
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
# rv key
RV_KEY = 'RV_OBJ'
OBJECTS = 'Gl699'

# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # -------------------------------------------------------------------------
    # ge bjd vs rv for all objects
    bjds = dict()
    rvs = dict()

    # loop around reductions
    for it, name in enumerate(paths):
        # get table
        table = Table.read(paths[name], format='ascii')
        # load into arrays
        bjds[name] = table['BJD']
        rvs[name] = table[RV_KEY] * 1000
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
