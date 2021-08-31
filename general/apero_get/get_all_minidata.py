#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2021-08-31

@author: cook
"""
from astropy.io import fits
import glob
import os
import shutil
from tqdm import tqdm


# =============================================================================
# Define variables
# =============================================================================
workdir = '/scratch2/drs-data/'
output = '/spirou2/minidata/'

targets = ['GL699']

# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":

    # get directories
    directories = glob.glob(workdir + 'minidata*')

    # loop around directories
    for directory in directories:
        # print progress
        print('Processing directory: {0}'.format(directory))
        # get outdir
        outdir = os.path.join(output, os.path.basename(directory))
        # rm outdir if it exists
        if os.path.exists(outdir):
            os.system('rm -rfv {0}'.format(outdir))
        # make out dir
        os.mkdir(outdir)
        # get pattern
        pattern = directory + '/reduced/'

        # get files
        e2ds_files = glob.glob(pattern + '*/2*o_pp_e2dsff_AB.fits')
        tcorr_files = glob.glob(pattern + '*/2*o_pp_e2dsff_tcorr_AB.fits')
        templates = glob.glob(pattern + 'other/Template_s1d_*AB.fits')

        all_files = e2ds_files + tcorr_files + templates

        # loop around files and add to outdir (if correct target)
        for filename in tqdm(all_files):
            # get header for this file
            header = fits.getheader(filename, ext=0)
            if 'DRSOBJN' in header:
                objname = str(header['DRSOBJN'])
            else:
                header = fits.getheader(filename, ext=1)
                objname = str(header['DRSOBJN'])
            del header
            # do not continue if objname not in targets
            if objname not in targets:
                continue
            # reate out path
            outpath = os.path.join(outdir, objname)
            # create object name dir
            if not os.path.exists(outpath):
                os.mkdir(outpath)
            # construct outpath
            outfilename = os.path.join(outpath, os.path.basename(filename))
            # copy
            shutil.copy(filename, outfilename)

# =============================================================================
# End of code
# =============================================================================
