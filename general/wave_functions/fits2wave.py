#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-06-22

@author: cook
"""
import numpy as np


# =============================================================================
# Define variables
# =============================================================================


# =============================================================================
# Define functions
# =============================================================================
def fits2wave(image, hdr):
    """
    Get the wave solution from the header using a filename
    """
    # size of the image
    nbypix, nbxpix = image.shape
    # get the keys with the wavelength polynomials
    wave_hdr = hdr['WAVE0*']
    # concatenate into a numpy array
    wave_poly = np.array([wave_hdr[i] for i in range(len(wave_hdr))])
    # get the number of orders
    nord = hdr['WAVEORDN']
    # get the per-order wavelength solution
    wave_poly = wave_poly.reshape(nord, len(wave_poly) // nord)
    # project polynomial coefficiels
    wavesol = np.zeros_like(image)
    # loop around orders
    for order_num in range(nord):
        wavesol[order_num] = np.polyval(wave_poly[order_num][::-1],np.arange(nbxpix))
    # return wave grid
    return wavesol




# =============================================================================
def fits2wavenew(image, hdr):
    """
    Get the wave solution from the header using a filename
    """
    # size of the image
    nbypix, nbxpix = image.shape
    # get the keys with the wavelength polynomials
    wave_hdr = hdr['WAVE0*']
    # concatenate into a numpy array
    wave_poly = np.array([wave_hdr[i] for i in range(len(wave_hdr))])
    # get the number of orders
    nord = hdr['WAVEORDN']
    # get the per-order wavelength solution
    wave_poly = wave_poly.reshape(nord, len(wave_poly) // nord)
    # project polynomial coefficiels
    wavesol = np.zeros_like(image)
    # xpixel grid
    xpix = np.arange(nbxpix)
    # loop around orders
    for order_num in range(nord):
        # calculate wave solution for this order
        owave = val_cheby(wave_poly[order_num], xpix, domain=[0, nbxpix])
        # push into wave map
        wavesol[order_num] = owave
    # return wave grid
    return wavesol


def val_cheby(coeffs, xvector,  domain):
    """
    Using the output of fit_cheby calculate the fit to x  (i.e. y(x))
    where y(x) = T0(x) + T1(x) + ... Tn(x)

    :param coeffs: output from fit_cheby
    :param xvector: x value for the y values with fit
    :param domain: domain to be transformed to -1 -- 1. This is important to
    keep the components orthogonal. For SPIRou orders, the default is 0--4088.
    You *must* use the same domain when getting values with fit_cheby
    :return: corresponding y values to the x inputs
    """
    # transform to a -1 to 1 domain
    domain_cheby = 2 * (xvector - domain[0]) / (domain[1] - domain[0]) - 1
    # fit values using the domain and coefficients
    yvector = np.polynomial.chebyshev.chebval(domain_cheby, coeffs)
    # return y vector
    return yvector


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # print hello world
    print('Hello World')

# =============================================================================
# End of code
# =============================================================================
