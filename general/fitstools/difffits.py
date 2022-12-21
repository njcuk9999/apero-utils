#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-12-19 at 10:29

@author: cook
"""
from astropy.io import fits
import sys
import os


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    _, file1, file2 = sys.argv

    data1 = fits.getdata(file1)
    data2 = fits.getdata(file2)

    difffile = '/home/cook/apero-data/diff.fits'

    fits.writeto(difffile, data1 - data2, overwrite=True)

    os.system(f'/scratch/bin/ds9/ds9 {file1} {difffile}')

# =============================================================================
# End of code
# =============================================================================
