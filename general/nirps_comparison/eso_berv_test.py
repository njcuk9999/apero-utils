#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-06-17 at 14:04

@author: cook
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from tqdm import tqdm
import warnings
import glob

# =============================================================================
# Define variables
# =============================================================================
YKEY = 'HIERARCH ESO QC TMMEAN USED'
XKEY = 'MJD-OBS'
TITLE = 'TMMEAN used in BERV'
# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def function1():
    return 0


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # print 'Hello World!'
    files = glob.glob('*/*_S2D_A.fits')

    tmmean = []
    mjd = []

    for filename in tqdm(files):
        hdr = fits.getheader(filename)
        tmmean.append(float(hdr[YKEY]))
        mjd.append(float(hdr[XKEY]))

    plt.close()
    fig, frames = plt.subplots(nrows=2, ncols=1)

    frames[0].plot(mjd, tmmean, marker='.', ls='None')
    frames[1].hist(tmmean, bins=10)

    frames[0].set(xlabel=XKEY, ylabel=YKEY)
    frames[1].set(xlabel=YKEY)

    plt.show()

# =============================================================================
# End of code
# =============================================================================
