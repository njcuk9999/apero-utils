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

# =============================================================================
# Define variables
# =============================================================================


# =============================================================================
# Define functions
# =============================================================================

files = glob.glob('qc_check/*f_pp.fits')
i1 = 0
im1, hdr = fits.getdata(files[i1], header=True)
print(hdr['UTC-OBS'],hdr['DPRTYPE'],files[i1])
im1 /= np.sqrt(np.nansum(im1 ** 2))
for i2 in range(1,len(files)):
    im2,hdr = fits.getdata(files[i2],header = True)
    im2 /= np.sqrt(np.nansum(im2**2))
    print(i1,i2,np.nansum(im1*im2),hdr['UTC-OBS'],hdr['DPRTYPE'],files[i2])


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    print('Hello World')

# =============================================================================
# End of code
# =============================================================================
