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
# INSTRUMENT = 'NIRPS_HA'
INSTRUMENT = 'NIRPS_HE'


if INSTRUMENT == 'NIRPS_HA':
    INDIR = '/nirps_raw/nirps/raw-data/nirps_ha/'
    OUTDIR = '/nirps_raw/nirps/apero-data/common/rawsym202205-HA/'
    NIGHTS = ['2022-05-17', '2022-05-18', '2022-05-19', '2022-05-20',
              '2022-05-21', '2022-05-22', '2022-05-23']
else:
    INDIR = '/nirps_raw/nirps/raw-data/nirps_he/'
    OUTDIR = '/nirps_raw/nirps/apero-data/common/rawsym202205-HE/'
    NIGHTS = ['2022-05-17', '2022-05-18', '2022-05-20',
              '2022-05-21', '2022-05-22', '2022-05-23']

DEBUG = False
SYMLINK = True


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":

    # get directories
    for night in NIGHTS:
        # construct indir
        indir = os.path.join(INDIR, night)
        # construct outdir
        outdir = os.path.join(OUTDIR, night)
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
            print('{0} --> {1}'.format(inpath, outpath))
            # if not in debug mode copy
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
