#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2021-07-08

@author: cook
"""
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import os
import sys


# =============================================================================
# Define variables
# =============================================================================
path = '/data/spirou/data/minidata2/reduced/2020-08-31/'



# =============================================================================
# Define functions
# =============================================================================
def diff_image(imagepath, imagename):
    try:
        hdu1 = fits.open(os.path.join(imagepath, imagename))
        hdu2 = fits.open(imagename)
    except:
        print('Skipping {0} [non-fits]'.format(imagename))
        return

    for extnum in range(len(hdu1)):

        # get name
        name = '{0}[{1}]'.format(imagename, extnum)
        print('=' * 50)
        print(name)
        print('=' * 50)

        if extnum >= len(hdu2):
            print('\tEXTENSION {0} MISSING HDU2'.format(extnum))
            continue

        # deal with image hdu
        if isinstance(hdu1[extnum], fits.ImageHDU):
            imdiff = fits.diff.ImageDataDiff(hdu1[extnum].data, hdu2[extnum].data)
            print(imdiff.report())

            diff = hdu1[extnum].data - hdu2[extnum].data
            if np.nansum(diff) != 0:
                fig, frame = plt.subplots(ncols=1, nrows=1)
                pos = frame.imshow(diff, aspect='auto', origin='lower')
                frame.set(title=name)
                fig.colorbar(pos, ax=frame)
                plt.show()
                plt.close()
        elif isinstance(hdu1[extnum], fits.BinTableHDU):
            imdiff = fits.diff.TableDataDiff(hdu1[extnum].data, hdu2[extnum].data)
            print(imdiff.report())
        else:
            print('\tSkipping (not ImageHDU or BinHDU)')


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    files = np.array(os.listdir('.'))
    last_modified = []
    # get last modified for all files
    for filename in files:
        last_modified.append(os.path.getmtime(filename))
    # sort by last modified
    sortmask = np.argsort(last_modified)
    files = files[sortmask]
    # diff images in order
    for filename in files:
        diff_image(path, filename)


# =============================================================================
# End of code
# =============================================================================
