#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-05-18

@author: cook
"""
from astropy.io import fits
import glob
import os
import shutil


# =============================================================================
# Define variables
# =============================================================================
INSTRUMENT = 'NIRPS_HA'
# INSTRUMENT = 'NIRPS_HE'


if INSTRUMENT == 'NIRPS_HA':
    INDIR = '/nirps_raw/nirps/raw-data/2022*'
    OUTDIR = '/nirps_raw/nirps/raw-data/nirps_ha/'
    MODE = 'HA'
else:
    INDIR = '/nirps_raw/nirps/raw-data/2022*'
    OUTDIR = '/nirps_raw/nirps/raw-data/nirps_he/'
    MODE = 'HE'

DEBUG = False

SYMLINK = False

FILTER = dict()
FILTER['HIERARCH ESO INS MODE'] = MODE


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

    # get directories
    directories = glob.glob(INDIR)


    for directory in directories:

        # construct indir
        indir = str(directory) + os.sep
        # construct outdir
        outdir = os.path.join(OUTDIR, os.path.basename(directory))

        # make output dir
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        # get file list
        files = os.listdir(indir)

        print('{0} files found'.format(len(files)))
        # loop around files
        for basename in files:
            # get paths
            inpath = os.path.join(indir, basename)
            outpath = os.path.join(outdir, basename)
            # remove inpath if outpath exists
            if os.path.exists(outpath):
                continue
            # skip non-fits files
            if not basename.endswith('.fits'):
                continue
            # get header
            hdr = fits.getheader(inpath)
            # filter files by header values in FILTER
            if filter_files(hdr):
                print('{0} --> {1}'.format(inpath, outpath))
                if not DEBUG:
                    try:
                        if SYMLINK:
                            os.symlink(inpath, outpath)
                        else:
                            shutil.move(inpath, outpath)
                    except Exception as _:
                        continue

# =============================================================================
# End of code
# =============================================================================
