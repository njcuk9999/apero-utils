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

# =============================================================================
# Define variables
# =============================================================================
# define raw directory
RAW_DIR = '/data/spirou/data/common/raw/'
# define out directory
OUT_DIR = '/data/spirou/data/common/rawfilesym/'


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # get all raw files
    files = glob.glob(os.path.join(RAW_DIR, '*', '*.fits'))
    # loop around files and make links
    for it, filename in enumerate(files):
        # get basename
        basename = os.path.basename(filename)
        # get directory
        obs_dir = os.path.basename(os.path.dirname(filename))
        # create outfile
        outfile = os.path.join(OUT_DIR, obs_dir, basename)
        # print progress
        margs = [it + 1, len(files), outfile]
        print('[{0}/{1}] {2}'.format(*margs))
        # make a symlink
        os.symlink(filename, outfile)


# =============================================================================
# End of code
# =============================================================================
