#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CODE DESCRIPTION HERE

Created on 2021-01-2021-01-08 15:29

@author: cook
"""
import numpy as np
import glob
from astropy.io import fits
import os
from scipy.stats import pearsonr

# =============================================================================
# Define variables
# =============================================================================


# =============================================================================
# Define functions
# =============================================================================

files = glob.glob('qc_check/*f_pp.fits')
i1 = 0
im1, hdr = fits.getdata(files[i1], header=True)
im1 = im1.ravel()
print(hdr['UTC-OBS'],hdr['DPRTYPE'],files[i1])

for i2 in range(1,len(files)):
    im2,hdr = fits.getdata(files[i2],header = True)
    im2 = im2.ravel()
    good = np.isfinite(im1) * np.isfinite(im2)
    metric, _ = pearsonr(im1[good], im2[good])
    print(i1,i2,metric,hdr['UTC-OBS'],hdr['DPRTYPE'],files[i2])


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    print('Hello World')

# =============================================================================
# End of code
# =============================================================================
