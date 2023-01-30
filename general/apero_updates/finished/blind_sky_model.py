#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-01-30 at 14:28

@author: cook
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from apero.core.math import lowpassfilter
from tqdm import tqdm
from apero.core.math import estimate_sigma as sigma
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from scipy.ndimage import binary_erosion


# =============================================================================
# Define variables
# =============================================================================

# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def function1():
    return 0


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------

    # retrieve A and B spectra with corresponding wavelengths
    wave_A = fits.getdata('NIRPS_2022-11-25T12_00_56_233_pp_e2dsff_A_wavesol_ref_A.fits')
    sp_A = np.array(fits.getdata('NIRPS_2022-12-06T01_30_31_491_pp_e2dsff_A.fits'), dtype=float)

    with fits.open('sky_model.fits') as hdulist:
        model_wave = hdulist['WAVE'].data
        model_A = hdulist['SCI_SKY'].data
        weights = hdulist['WEIGHTS'].data
        reg_id = hdulist['REG_ID'].data
        grad_model = hdulist['GRADIENT'].data

    # placeholder for the sky reconstructed for this file
    sky_A = np.zeros_like(sp_A)

    # TODO use proper wave shift
    # shift to the reference grid
    sp_A[~np.isfinite(sp_A)] = 0.0
    for iord in range(sp_A.shape[0]):
        valid = np.isfinite(sp_A[iord])
        sp_A[iord] = ius(wave_A[iord][valid], sp_A[iord][valid], ext=1, k=3)(model_wave[iord])

    # TODO subtract template here. Each order is scaled to its median. OH lines
    # TODO contribute very little to the median flux, so this is fine.

    # find gradients in file
    grad = np.gradient(sp_A, axis=1)

    for i in tqdm(np.unique(reg_id)[1:], leave=False):
        # find pixels that are within the region ID and valid
        g = (reg_id == i) * np.isfinite(model_A) * np.isfinite(sp_A) * \
            np.isfinite(grad) * np.isfinite(grad_model)

    grad1 = grad[g]
    grad2 = grad_model[g]
    # dot-product of derivatives
    amp = np.nansum(grad1 * grad2) / np.nansum(grad2 ** 2)
    sky_A[g] = model_A[g] * amp * weights[g]

    # TODO use proper wave shift
    # spline back to the file's wavelength grid
    for iord in range(sp_A.shape[0]):
        sky_A[iord] = ius(model_wave[iord], sky_A[iord], ext=1, k=3)(wave_A[iord])

    # some plots for fun
    fig, ax = plt.subplots(nrows=2, ncols=1, sharex='all', figsize=[16, 8])
    ax[0].plot(wave_A.ravel(), sp_A.ravel(), alpha=0.5)
    ax[0].plot(wave_A.ravel(), sp_A.ravel() - sky_A.ravel())
    ax[1].plot(wave_A.ravel(), sp_A.ravel() - sky_A.ravel())
    # ax.set(xlim =[1500,1700],ylim = [0,5000])
    ax[0].set(xlabel='Wavelength [nm]', ylabel='Normalized flux',
              title='Blind sky subtraction')
    plt.tight_layout()
    plt.savefig('sky_blind.png')
    plt.show()

    fits.writeto('sky_subtract_A.fits', sp_A - sky_A, overwrite=True)


# =============================================================================
# End of code
# =============================================================================
