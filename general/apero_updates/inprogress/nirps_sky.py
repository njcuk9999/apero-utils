#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-01-30 at 14:33

@author: cook
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from apero.core.math import lowpassfilter
import glob
from tqdm import tqdm
from apero.core.math import estimate_sigma
from scipy.ndimage import binary_erosion, binary_dilation
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from astropy.table import Table
from scipy.optimize import curve_fit
from scipy import constants

# =============================================================================
# Define variables
# =============================================================================

# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def supergauss(x, center, amp, fwhm, expo, zp):
    """
    :param x:
    :param center:
    :param amp:
    :param fwhm:
    :param expo:
    :param zp:
    :return:
    0.5 = exp(-0.5*((fwhm/2)/ew)^expo)
    np.log(0.5) = -0.5*((fwhm/2)/ew)^expo
    2*np.log(2) = ((fwhm/2)/ew)^expo
    """
    ew = (fwhm / 2) / (2 * np.log(2)) ** (1 / expo)
    return np.exp(-0.5 * np.abs((x - center) / ew) ** np.abs(expo)) * amp + zp


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # This code loads and processes sky data from NIRPS to identify and mask out absorption lines.
    # This is useful for removing sky lines that could affect the accuracy of spectral analysis.

    
    fiber = 'HE'
    doplot = True
    if fiber == 'HE':
        # get wavelength grid for A spectra
        waveref = fits.getdata('NIRPS_2022-11-25T12_00_56_233_pp_e2dsff_A_wavesol_ref_A.fits')
    
    # get reference A files
    all_sky_A = np.array(glob.glob('sky_HE/NIRPS_202*_pp_e2dsff_A.fits'))
    
    if fiber == 'HA':
        # get wavelength grid for A spectra
        waveref = fits.getdata('NIRPS_2022-11-25T12_24_04_050_pp_e2dsff_A_wavesol_ref_A.fits')
    
    # get reference A files
    all_sky_A = np.array(glob.glob('sky_HA/*_pp_e2dsff_A.fits'))
    
    # get corresponding reference B files
    all_sky_B = np.array([f.replace('_A.', '_B.') for f in all_sky_A])
    
    # find number of orders
    nord = fits.getdata(all_sky_A[0]).shape[0]
    
    Nbin = 5
    
    # place holder for the A and B fiber sky 'cubes'
    all_sky_map_A = np.zeros([Nbin, len(waveref) * 4088])
    all_sky_map_B = np.zeros_like(all_sky_map_A)
    
    index = (Nbin * np.arange(len(all_sky_A))) // len(all_sky_A)
    
    # loading all sky data
    for i in tqdm(range(len(all_sky_A))):
        sky_A = np.array(fits.getdata(all_sky_A[i]), dtype=float)
        hdr_A = fits.getheader(all_sky_A[i])
        wave_A = hdr2wave(hdr_A)
        # avoid NaNs for the spline later
        sky_A[~np.isfinite(sky_A)] = 0.0
    
        sky_B = np.array(fits.getdata(all_sky_B[i]), dtype=float)
        hdr_B = fits.getheader(all_sky_B[i])
        wave_B = hdr2wave(hdr_B)

        # avoid NaNs for the spline later
        sky_B[~np.isfinite(sky_B)] = 0.0

        for iord in range(nord):
            # high-pass the sky
            sky_A[iord] -= lowpassfilter(sky_A[iord], 101)
            sky_B[iord] -= lowpassfilter(sky_B[iord], 101)

        # spline onto the master grid
        sky_A[iord] = ius(wave_A[iord], sky_A[iord], ext=1, k=3)(waveref[iord])
        sky_B[iord] = ius(wave_B[iord], sky_B[iord], ext=1, k=3)(waveref[iord])

        # expressed as a ravelled e2ds
        all_sky_map_A[index[i]] += sky_A.ravel()
        all_sky_map_B[index[i]] += sky_B.ravel()
    
    # get ravel-e2ds of the wavelegnth grid for dilation stuff later
    waveref2 = waveref.ravel()
    
    # get median sky spectrum, only used for identification of lines
    v = np.nanmedian(all_sky_map_A, axis=0)
    v[~np.isfinite(v)] = 0
    
    # find positive excursions in sky signal
    nsig = v / sigma(v)
    # lines are >5 sigma positive excursion
    line = np.array(nsig > 5, dtype=int)
    # erode features that are too narrow
    line = binary_erosion(line, structure=np.ones(5))
    # dilate to get wings of lines
    line = binary_dilation(line, structure=np.ones(27))
    
    reg_id = (np.cumsum(line != np.roll(line, 1)))
    reg_id[(reg_id % 2) == 0] = 0
    reg_id[reg_id != 0] = (reg_id[reg_id != 0] + 1) // 2
    
    # put the line mask onto a magic grid to avoid errors at order overlaps
    # TODO use the magic grid function
    magic_grid = np.array(Table.read('NIRPS_2022-12-02T00_59_50_241_pp_s1d_v_A.fits')['wavelength'])
    magic_mask = np.zeros_like(magic_grid, dtype=bool)
    
    # put everything in the magic grid referential
    for i in np.unique(reg_id)[1:]:
        g = (i == reg_id)
        magic_mask[(magic_grid > np.min(waveref2[g])) * (magic_grid < np.max(waveref2[g]))] = True
    
    # now in the space of the magic grid
    reg_id_magic = (np.cumsum(magic_mask != np.roll(magic_mask, 1)))
    reg_id_magic[(reg_id_magic % 2) == 0] = 0
    reg_id_magic[reg_id_magic != 0] = (reg_id_magic[reg_id_magic != 0] + 1) // 2
    
    # fill the map with unique values and common ID for overlapping orders
    reg_id = np.zeros_like(reg_id, dtype=int)
    for ureg in np.unique(reg_id_magic)[1:]:
        wave_min = np.min(magic_grid[reg_id_magic == ureg])
        wave_max = np.max(magic_grid[reg_id_magic == ureg])
    
    g = (waveref2 > wave_min) * (waveref2 < wave_max)
    reg_id[g] = ureg
    
    
    # plots that you can remove, for debug
    
    vv = v.reshape(waveref.shape)
    for iord in range(waveref.shape[0]):
        plt.plot(waveref[iord], vv[iord], alpha=0.5)
    plt.plot(waveref2[reg_id != 0], v[reg_id != 0], 'g.', alpha=0.3)
    
    for reg in np.unique(reg_id):
        if reg == 0:
            continue
        gg = reg_id == reg
        plt.text(np.mean(waveref2[gg]), 0, '{}'.format(reg))
    plt.show()
    
    fig, ax = plt.subplots(nrows=2, ncols=1)
    # construct a model with all lines normalized to a median of 1 in fiber A
    model_A = np.zeros_like(v)
    model_B = np.zeros_like(v)
    
    fwhm_all = []
    xpix_all = []
    nsig_all = []
    iord, xpix = np.indices(waveref.shape)
    xpix = xpix.ravel()
    iord = iord.ravel()
    for ii in tqdm(np.unique(reg_id)[1:]):
        gg = reg_id.ravel() == ii
    
        all_A = []
        all_B = []
        for i in range(all_sky_map_B.shape[0]):
            tmp_A = all_sky_map_A[i][gg]
            tmp_B = all_sky_map_B[i][gg]

            amp = np.nansum(tmp_A)
            all_A = np.append(all_A, tmp_A / amp)
            all_B = np.append(all_B, tmp_B / amp)
    
            if doplot:
                ax[0].plot(waveref2[gg], tmp_A / amp, alpha=0.5, color='orange')
                ax[1].plot(waveref2[gg], tmp_B / amp, alpha=0.5, color='orange')

        med_A = np.nanmedian(all_A.reshape(Nbin, len(tmp_A)), axis=0)
        med_B = np.nanmedian(all_B.reshape(Nbin, len(tmp_B)), axis=0)
    
        try:
            dv = (waveref2[gg] / np.mean(waveref2[gg]) - 1) * (constants.c / 1e3)
            # get the FWHM of the line
            p0 = [0, np.max(med_A), 4.0, 2.0, 0.0]
            fit, cov = curve_fit(supergauss, dv, med_A, p0=p0)
            fwhm_all = np.append(fwhm_all, fit[2])
            xpix_all = np.append(xpix_all, np.mean(xpix[gg]))
            nsig = fit[1] / np.nanstd(med_A - supergauss(dv, *fit))
            nsig_all = np.append(nsig_all, nsig)
        except:
            _ = True

        if ii == 145:
            doplot = True
        else:
            doplot = False
        if doplot:
            fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True)


        if doplot:
            ax[0].plot(waveref2[gg], med_A, color='blue')
            ax[1].plot(waveref2[gg], med_B, color='blue')

        model_A[gg] = med_A
        model_B[gg] = med_B
        if doplot:
            plt.show()

    keep = nsig_all > 5
    xpix_all = xpix_all[keep]
    fwhm_all = fwhm_all[keep]
    nsig_all = nsig_all[keep]
    
    plt.plot(xpix_all, fwhm_all, 'go', alpha=0.5)
    plt.show()
    
    model_A = model_A.reshape(sky_A.shape)
    model_B = model_B.reshape(sky_B.shape)
    reg_id = reg_id.reshape(sky_A.shape)
    
    
    
    weights = np.array(reg_id != 0, dtype=float).ravel()
    weights2 = np.array(reg_id != 0, dtype=float).ravel()
    
    for i in range(5):
        weights2 = binary_erosion(weights2, structure=np.ones(3))
        print(np.mean(weights2))
        weights += np.array(weights2, dtype=float)
    
    weights /= np.max(weights)
    weights = weights.reshape(model_A.shape)
    
    hdr = fits.Header()
    hdr['HAS_SKY'] = True  # that's for NIRPS, SPIRou would have False
    hdu0 = fits.PrimaryHDU()
    hdulist = [hdu0]
    
    # set to 0 negative values that are spuriuous for a sky flux
    model_A[model_A < 0] = 0.0
    model_B[model_B < 0] = 0.0
    
    hdr = fits.Header()
    hdr['EXTNAME'] = 'SCI_SKY'
    hdulist.append(fits.ImageHDU(data=model_A, header=hdr))
    
    hdr = fits.Header()
    hdr['EXTNAME'] = 'CAL_SKY'
    # TODO fill with zeros for SPIRou and SPIP and VROOMM
    hdulist.append(fits.ImageHDU(data=model_B, header=hdr))
    
    hdr = fits.Header()
    hdr['EXTNAME'] = 'WAVE'
    hdulist.append(fits.ImageHDU(data=waveref, header=hdr))
    
    hdr = fits.Header()
    hdr['EXTNAME'] = 'REG_ID'
    hdulist.append(fits.ImageHDU(data=reg_id, header=hdr))
    
    hdr = fits.Header()
    hdr['EXTNAME'] = 'WEIGHTS'
    hdulist.append(fits.ImageHDU(data=weights, header=hdr))
    
    hdr = fits.Header()
    hdr['EXTNAME'] = 'GRADIENT'
    hdulist.append(fits.ImageHDU(data=np.gradient(model_A, axis=1), header=hdr))
    
    hdul = fits.HDUList(hdulist)
    hdul.writeto('sky_model_{}.fits'.format(fiber.lower()), overwrite=True)
# =============================================================================
# End of code
# =============================================================================
