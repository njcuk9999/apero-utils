#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2021-06-14

@author: cook
"""
import os
import glob
from tqdm import tqdm

# =============================================================================
# Define variables
# =============================================================================
# define raw directory
RAW_DIR = '/spirou/cfht_nights/common/raw'
# define out directory
OUT_DIR = '/spirou/drs-data/common/raw'


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # get all raw files
    files = glob.glob(os.path.join(RAW_DIR, '*', '*.fits'))
    # loop around files and make links
    for it, filename in tqdm(enumerate(files)):
        # get basename
        basename = os.path.basename(filename)
        # get directory
        obs_dir = os.path.basename(os.path.dirname(filename))
        # create outfile
        outfile = os.path.join(OUT_DIR, obs_dir, basename)
        # do not recreate links
        if not os.path.exists(outfile):
            # create dir if it doesn't exist
            if not os.path.exists(os.path.dirname(outfile)):
                os.mkdir(os.path.dirname(outfile))
            # print progress
            margs = [it + 1, len(files), outfile]
            print('[{0}/{1}] {2}'.format(*margs))
            # make a symlink
            os.symlink(filename, outfile)


# =============================================================================
# End of code
# =============================================================================
