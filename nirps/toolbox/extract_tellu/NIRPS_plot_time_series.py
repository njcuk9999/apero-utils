#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to plot a detrended spectral time series, akin to what is done for WP3 analyses. 
Input: list is NIRPS files
Output: plots
"""
import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits


# list containing the path to all NIRPS exposures to load
file_list = []


#%% DEFINE FUNCTIONS

# functions from Neil to read NIRPS wavelength from header keys
def val_cheby(coeffs, xvector, domain):
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

def read_NIRPS_wave(header):
    """
        Provide a fits header
        and get the corresponding wavelength
        grid from the header.
        Usage :
          wave = fits2wave(hdr)
                  or
          wave = fits2wave('my_e2ds.fits')
        Output has the same size as the input
        grid. This is derived from NAXIS
        values in the header
    """
    # get the keys with the wavelength polynomials
    wave_hdr = header['WAVE0*']
    # concatenate into a numpy array
    wave_poly = np.array([wave_hdr[i] for i in range(len(wave_hdr))])
    # get the number of orders
    try:
        nord = header['WAVEORDN']
    except:
        nord = header['NAXIS2']
    # get the per-order wavelength solution
    wave_poly = wave_poly.reshape(nord, len(wave_poly) // nord)
    # get the length of each order (normally that's 4088 pix)
    try:
        npix = header['NAXIS1']
    # APERO 0.7 moved 'NAXIS1' to second header
    except:
        npix = 4088
    # project polynomial coefficiels
    wavesol = [val_cheby(wave_poly[i],np.arange(npix), domain=[0, npix]) for i in range(nord) ]
    # return wave grid
    return np.array(wavesol)

#%% LOAD DATA

# define quantities to be extracted as a list
counts_start = []

# loop over each exposure in file list
print('Loading Exposures')
for ii, exposure in enumerate(file_list):
    # open file
    hdul = fits.open(data_path + exposure)
    # read header
    file_header = hdul[0].header
    
    # extract flux
    exposure_counts = hdul[1].data

    # append values to keep in lists
    counts_start.append(exposure_counts)

# read wavelength array from header keys 
wave = read_NIRPS_wave(file_header) /1e3 # change to microns

nexposures = len(file_list)
norders, npixels = wave.shape
orders = np.arange(norders)

# convert all lists into numpy arrays
counts_start = np.array(counts_start)


#%% PROCESS DATA

# list that will keep track of intermediate steps in the detrending
reduction_steps = []
step_labels = []

# copy start data
flux = counts_start*1

# start data after bad columns removed
reduction_steps.append(flux*1)
step_labels.append('Start data')

# set all columns with any value below 0.5 to nans
flux[flux < 0.5] = np.nan
means = np.mean(flux,axis=0)
for order in orders:
    flux[:,order,np.isnan(means[order])] = np.nan

# normalize
flux = flux/np.nanmedian(flux, axis=(1,2))[:,None,None]
step_labels.append('After continuum normalization')

# data after normalization
reduction_steps.append(flux*1)

flux /= np.nanmedian(flux, axis=0)

# data after median spectrum removal
reduction_steps.append(flux*1)
step_labels.append('After median spectrum removed')




#%% PLOT DATA

# choose spectral order to plot
order = 37

XX1, YY1 = np.meshgrid(np.arange(nexposures), wave[order,:])

fig, ax = plt.subplots(len(reduction_steps), 1, sharex=True, figsize=(8,6))

for i in range(len(reduction_steps)):
            
    ax[i].tick_params(labelbottom=None)
    ax[i].set_title(step_labels[i])
            
    im = ax[i].pcolormesh(YY1,XX1, np.transpose(reduction_steps[i][:,order,:]), cmap='viridis', vmin = np.nanpercentile(reduction_steps[i][:,order,:],1), vmax= np.nanpercentile(reduction_steps[i][:,order,:],99),rasterized=True )
    ax[i].set_rasterization_zorder(-10)
    cb = plt.colorbar(im, ax=ax[i], pad=0.01, aspect=10)
    cb.minorticks_off()

    ax[i].set_ylabel('Exposure #')

    if i == int(len(reduction_steps)-1):
        ax[i].tick_params(labelbottom=True, labeltop=None, right = True)

    ax[i].minorticks_off()
    
ax[-1].set_xlabel(r'Wavelength [$\mu$m]')


fig.tight_layout()
# Remove horizontal space between axes
fig.subplots_adjust(hspace=0.2)










