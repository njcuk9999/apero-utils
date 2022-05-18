#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-05-18

@author: cook
"""
from astropy.io import fits
import os

# =============================================================================
# Define variables
# =============================================================================
indir = '/nirps_raw/nirps/raw-data/2022-05-17'
outdir = '/nirps_raw/nirps/apero-data/common/rawsym202205-HA/2022-05-17'

DEBUG = False

FILTER = dict()
FILTER['HIERARCH ESO INS MODE'] = 'HA'

# =============================================================================
# Define functions
# =============================================================================
def filter_files(hdr):
    """
    Filter files by header keys
    """
    for key in FILTER:
        if key in hdr:
            if hdr[key] != FILTER[key]:
                return False
    return True

# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # make output dir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # get file list
    files = os.listdir(indir)
    # loop around files
    for basename in files:
        # get paths
        inpath = os.path.join(indir, basename)
        outpath = os.path.join(outdir, basename)
        # skip non-fits files
        if not basename.endswith('.fits'):
            continue
        # get header
        hdr = fits.getheader(inpath)
        # filter files by header values in FILTER
        if filter_files(hdr):
            print('{0} --> {1}'.format(inpath, outpath))
            if not DEBUG:
                os.symlink(inpath, outpath)

# =============================================================================
# End of code
# =============================================================================
