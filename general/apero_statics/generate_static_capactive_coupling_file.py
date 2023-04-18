#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-04-18 at 13:14

@author: cook
"""

from astropy.io import fits
from etienne_tools import robust_polyfit,lowpassfilter
from tqdm import tqdm

import numpy as np
import matplotlib.pyplot as plt
import glob
import os

model_name = 'amplifier_bias_model_nirps.fits'

if not os.path.isfile(model_name):
    # find all dark files and construct ribbons for all of them
    files = np.array(glob.glob('NIRPS_PATCH*.fits'))

    for file in files:
        # if the ribbon for that file exists, we skip
        outname = 'ribbon_'+file
        if not os.path.isfile(outname):
            # read the dark frame
            dark = fits.getdata(file)
            # read its header
            h = fits.getheader(file)
            # subtract the per-column median to remove bulk DC level
            for icol in tqdm(range(dark.shape[0]),leave = False):
                dark[:,icol] -= np.nanmedian(dark[:,icol])


            namp = 32 # number of amplifier
            # pixel width of each amplifier
            pix_amp = dark.shape[1]//namp
            # cube of all amplifiers
            cube = np.zeros([dark.shape[0],pix_amp,namp])
            for ibin in tqdm(range(namp),leave = False):
                #  logics for the butterfly symmetry of amplifiers
                if (ibin % 2) == 0:
                    flip = 1
                else:
                    flip = -1
                # fill the amplifier cube with amps in the same directoin
                cube[:,:,ibin] = dark[:,pix_amp*ibin:pix_amp*(ibin+1)][:,::flip]

            # get the median amplifier pattern for file and save
            ribbon_pattern = np.nanmedian(cube,axis=2)
            fits.writeto(outname,ribbon_pattern,h)

    # find all ribbons
    files = np.array(glob.glob('amplifier_bias_model*.fits'))
    # get size of each ribbon
    sz = fits.getdata(files[0]).shape
    # build a big cube of these and keep track of exposure times
    cube = np.zeros([sz[0],sz[1],len(files)])
    exptime = np.zeros(len(files))
    for ifile in range(len(files)):
        h = fits.getheader(files[ifile])
        cube[:,:,ifile] = fits.getdata(files[ifile])
        exptime[ifile] = h['EXPTIME']

    # we expect the structure to scale with accumulated photons but express
    # everything as slopes. We thefore fit values linearly with
    # pix value * exptime
    slope = np.zeros(sz)
    intercept = np.zeros(sz)
    # for each ribbon pixel, we find the slope and intercept of the dark
    print('We fit intercept and slope')
    for i in tqdm(range(sz[0])):
        for j in range(sz[1]):
            fit,_ = robust_polyfit(exptime,cube[i,j,:]*exptime,1,3)
            slope[i,j] = fit[0]
            intercept[i,j] = fit[1]
    # save into a 2-extension fits file
    hdu0 = fits.PrimaryHDU()
    hdr = fits.Header()
    hdr['EXTNAME'] = 'slope'
    hdu1 = fits.ImageHDU(slope,header = hdr)

    hdr = fits.Header()
    hdr['EXTNAME'] = 'intercept'
    hdu2 = fits.ImageHDU(intercept,header = hdr)

    hdu = fits.HDUList([hdu0,hdu1,hdu2])
    hdu.writeto(model_name,overwrite = True)

# -----------------------------------------------------------------------------
# AFTER HERE IS JUST A TEST TO SEE THAT IT WORKS
# read a random file
file = 'NIRPS_PATCHED_2022-12-01T19_11_01_158.fits'
data = fits.getdata(file)
h = fits.getheader(file)

# read slope and intercept of model
slope = fits.getdata(model_name,extname = 'slope')
intercept = fits.getdata(model_name,extname = 'intercept')

# construct model for given file
amp_model = slope+intercept/h['EXPTIME']

# unfold the butterfly pattern
corr = np.zeros_like(data)
namp = 32
pix_amp = data.shape[1] // namp
for ibin in range(namp):
    if (ibin % 2) == 0:
        flip = 1
    else:
        flip = -1
    corr[:,pix_amp * ibin:pix_amp * (ibin+1)] =  amp_model[:,::flip]

# apply correction
data -= corr

fits.writeto('test_2d.fits',corr,h,overwrite = True)
fits.writeto('diff_2d.fits',data,h,overwrite = True)

# =============================================================================
# End of code
# =============================================================================
