#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-01-16 at 11:12

@author: cook
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from tqdm import tqdm
import warnings


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
    from astropy.io import fits
    from scipy.ndimage import binary_dilation
    import numpy as np
    from etienne_tools import sigma
    from etienne_tools import lin_mini
    from scipy.ndimage import median_filter
    import matplotlib.pyplot as plt


    def get_butterfly_maps(data, amp_flux=1.0, amp_deriv=1.0, amp_deriv2=1.0,
                           fast=False):
        """
        Usage --
        :param data: Raw NIRPS or SPIRou image for which we want to derive the
                     butterfly pattern arrising from capacitive compling. The
                     image should be an output of the fits2ramp and *not* a _pp file

        :param amp_flux: amplitude of the flux-dependend component

        :param amp_deriv: amplitude of the flux-dependent along-readout-axis
                          derivative component

        :param amp_deriv2: amplitude of the flux-dependent along-readout-axis
                          2nd derivative component
        :param fast: If true, only compute the total map

        :return: 4 images ---
                    - the full butterfly pattern
                    - the amplitude butterfly pattern [None if fast = True]
                    - the 2 derivatives butterfly pattern [None if fast = True]

        This function can be used in two ways; first to fit the find the default
        amplitude for a given science array (use the amplitude and derivative
        patterns) or, if you have these amplitudes from a previous file, as a
        predictor of the pattern for all science files.
        """
        namp = 32  # number of amplifier
        # pixel width of each amplifier
        pix_amp = data.shape[1] // namp
        # cube of all amplifiers
        cube = np.zeros([data.shape[0], pix_amp, namp])
        for ibin in range(namp):
            #  logics for the butterfly symmetry of amplifiers
            if (ibin % 2) == 0:
                flip = 1
            else:
                flip = -1
            # fill the amplifier cube with amps in the same directoin
            cube[:, :, ibin] = data[:, pix_amp * ibin:pix_amp * (ibin + 1)][:, ::flip]

        avg_amp = -np.nansum(cube, axis=2)
        _, grady = np.gradient(avg_amp)
        _, grady2 = np.gradient(grady)

        # we do not compute the full image of each component and do the sum on the
        # ribbons
        if fast:
            map_fast = amp_flux * avg_amp + amp_deriv * grady + amp_deriv2 * grady2
            full_butterfly = np.zeros_like(data)

            map_fast_flip = map_fast[:, ::-1]

            flux_buttefly = None
            deriv_butterfly = None
            deriv2_butterfly = None
        else:
            flux_buttefly = np.zeros_like(data)
            deriv_butterfly = np.zeros_like(data)
            deriv2_butterfly = np.zeros_like(data)

            avg_amp_flip = avg_amp[:, ::-1]
            grady_flip = grady[:, ::-1]
            grady2_flip = grady2[:, ::-1]

        for ibin in range(namp):
            # pixels defining the region of the amplifier
            i1 = pix_amp * ibin
            i2 = pix_amp * (ibin + 1)
            #  logics for the butterfly symmetry of amplifiers
            if not fast:
                if (ibin % 2) == 0:
                    flux_buttefly[:, i1:i2] = avg_amp
                    deriv_butterfly[:, i1:i2] = grady
                    deriv2_butterfly[:, i1:i2] = grady2
                else:
                    flux_buttefly[:, i1:i2] = avg_amp_flip
                    deriv_butterfly[:, i1:i2] = grady_flip
                    deriv2_butterfly[:, i1:i2] = grady2_flip
            else:
                if (ibin % 2) == 0:
                    full_butterfly[:, i1:i2] = map_fast
                else:
                    full_butterfly[:, i1:i2] = map_fast_flip

        if not fast:
            # construct the full butterfly pattend
            full_butterfly = amp_flux * flux_buttefly + \
                             amp_deriv * deriv_butterfly + amp_deriv2 * deriv2_butterfly

        return full_butterfly, flux_buttefly, deriv_butterfly, deriv2_butterfly


    def get_butterfly_amplitudes(data0, n_ite=6, sigma_cut=2.0):
        """

        :param data0: A raw image used to derive the amplitude of the flux and
                      derivative components of the butterfly pattern. This should be
                      a DARK+FLAT or DARK+FLAT for either SPIRou or NIRPS. Use a
                      file with a lot of flux but not saturated
                      For NIRPS do not use HE use HA

        :param flux_buttefly: The flux_butterfly output from 'get_butterfly_maps'
        :param deriv_butterfly: The deriv_butterfly output from 'get_butterfly_maps'
        :param deriv2_butterfly: The 2nd deriv_butterfly output
        :param n_ite: The number of sigma-clipping iterations. Leave default at 6
                      unless you have a good reason.
        :param sigma_cut: sigma-clipping threshold. Leave at 2.0 unless you have a
                          good reason to do so.
        :return: best-fitting (in a least-square sense) amplitudes from the flux
                 and derivative components of butterfly pattern
        """

        # we'll pad data with NaNs, to avoid shallow-copy errors, we create another
        # array
        data = np.array(data0, dtype=float)

        full_butterfly, flux_buttefly, deriv_butterfly, deriv2_butterfly = \
            get_butterfly_maps(data)

        # we will first create a version of the image where the out-of-order domain
        # is close to zero. This allows us to create a mask to find unillunimated
        # pixels.
        for ite in range(3):
            # iteratively subtraction median of columns and rows
            med = np.nanmedian(data, axis=0)
            med = np.tile(med, data.shape[0]).reshape(data.shape)
            data -= med
            med = np.nanmedian(data, axis=1)
            med = np.repeat(med, data.shape[0]).reshape(data.shape)
            data -= med
        # create a mask of points at >1 sigma. These are illuminated pixels
        mask = data > sigma(data)
        # dilate region to catch diffused light
        mask = binary_dilation(mask, np.ones([7, 7]))

        data[mask] = np.nan
        # as we are measuring a small signal correlated in the dispersion direction
        # we filter the image to increase the contribution of these spatial
        # structures are their measurement accuracy
        data = median_filter(data, [7, 1])

        valid = ~mask
        for ite in range(n_ite):
            # we construct a linear system and add a constant DC value, which does
            # not count in reducing dispersion but accounts for possible DC offset
            # of the image
            sample = np.array([flux_buttefly[valid], deriv_butterfly[valid],
                               deriv2_butterfly[valid], np.ones_like(data[valid])])
            fit, recon = lin_mini(data[valid], sample)

            diff = data[valid] - recon
            sig = sigma(diff)
            nsig = np.zeros_like(data) + np.inf
            nsig[valid] = np.abs(diff) / sig

            txtout = 'iteration {}/{}, improvement in rms {:5.2f}%, ' \
                     'flux amp {:2.4e}, deriv amp {:2.4e}, deriv2 amp {:2.4e}, ' \
                     'DC offset {:2.2e}'
            txtout = txtout.format(ite + 1, n_ite, (1 - sig / sigma(data[valid])) * 100,
                                   fit[0], fit[1], fit[2], fit[3])
            print(txtout)

            valid = np.isfinite(data) * (nsig < sigma_cut)

        print('\n\t ' + 60 * '*' + '\n')
        print('\t \t\t Recommended APERO default values\n')
        print('\t \t amp_flux = {:2.6e} '.format(fit[0]))
        print('\t \t amp_deriv = {:2.6e} '.format(fit[1]))
        print('\t \t amp_deriv2 = {:2.6e} '.format(fit[2]))
        print('\n\t ' + 60 * '*' + '\n')

        return fit[0], fit[1], fit[2]


    # sample NIRPS file
    # sample FLAT,DARK,LAMP
    # define your coefficients with this one
    # file = 'NIRPS_2022-11-25T11_46_49_037.fits'

    # sample GJ3090 file to apply coefficients
    # file = 'NIRPS_2022-12-02T03_56_36_492.fits'

    # sample SPIRou file
    # sample DARK,FLAT
    # define your coefficients with this one
    file = '2691856f.fits'

    data0 = fits.getdata(file)
    data = np.array(data0)

    h = fits.getheader(file)

    if True:
        # demo to measure the amplitudes. Use only with a DARK+FLAT (SPIRou)
        # or either DARK+FLAT or FLAT+DARK (HA) for NIRPS
        amp_flux, amp_deriv, amp_deriv2 = \
            get_butterfly_amplitudes(data0, n_ite=10)

    # Default NIRPS values
    amp_flux = 1.359371e-04
    amp_deriv = 7.727465e-04
    amp_deriv2 = -2.601081e-04

    # Default SPIRou values
    amp_flux = 1.490711e-04
    amp_deriv = 8.080468e-04
    amp_deriv2 = -6.523985e-04

    full_butterfly, _, _, _ = get_butterfly_maps(data, amp_flux, amp_deriv,
                                                 amp_deriv2, fast=True)

    fits.writeto('cleaned_{}'.format(file), data0 - full_butterfly, overwrite=True)

# =============================================================================
# End of code
# =============================================================================
