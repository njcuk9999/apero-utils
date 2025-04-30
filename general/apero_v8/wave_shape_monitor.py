#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-04-16 at 13:46

@author: cook
"""
import numpy as np
from astropy.io import fits
from tqdm import tqdm
import os
import glob
from astropy import constants
import matplotlib.pyplot as plt
from astropy.time import Time


# =============================================================================
# Define variables
# =============================================================================
PATH = '/cosmos99/nirps/apero-data/nirps_ha_online/red'
PATH = '/cosmos99/nirps/apero-data/nirps_he_online/red'
FILENAME = '*/*wave_night_A.fits'

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
    # get the wave files
    wavefiles = glob.glob(os.path.join(PATH, FILENAME))

    mjdmid = []
    wave_cents = []
    # loop around
    for wavefile in tqdm(wavefiles):
        # load the file
        wavemap = fits.getdata(wavefile)
        wavehdr = fits.getheader(wavefile)
        # get central pixel
        cent_x = wavemap.shape[1] // 2
        # get the mjd mid from the header
        mjdmid.append(wavehdr['MJDMID'])
        # get the central pixels of every order
        wave_cents.append(wavemap[:, cent_x])


    mjdmid = Time(np.array(mjdmid), format='mjd')
    wave_cents = np.array(wave_cents)

    dv = np.log(wave_cents / np.nanmedian(wave_cents, axis=0)) * constants.c.value

    dv_15 = np.nanmean(dv[:, 10:20], axis=1)
    dv_60 = np.nanmean(dv[:, 55:65], axis=1)

    mad_15 = np.nanmedian(abs(np.diff(dv_15))) / np.sqrt(2)
    mad_60 = np.nanmedian(abs(np.diff(dv_60))) / np.sqrt(2)


    fig, frames = plt.subplots(nrows=2, ncols=1, sharex='all')


    frames[0].plot_date(mjdmid.plot_date, dv_15, label='Order 15', ls='None', marker='o',
             color='purple')
    frames[0].plot_date(mjdmid.plot_date, dv_60, label='Order 60', ls='None', marker='x',
             color='orange')

    frames[0].set(xlabel='Date', ylabel='Central pixel offset [m/s]')

    frames[0].legend(loc=0)


    frames[1].plot_date(mjdmid.plot_date[1:], np.diff(dv_15) / np.diff(mjdmid.mjd), label='Order 15', ls='None', marker='o',
             color='purple')
    frames[1].plot_date(mjdmid.plot_date[1:], np.diff(dv_60) / np.diff(mjdmid.mjd), label='Order 60', ls='None', marker='x',
             color='orange')

    frames[1].set(xlabel='Date', ylabel='Central pixel offset [m/s/day]')

    frames[1].legend(loc=0)
    plt.show()
    plt.close()

# =============================================================================
# End of code
# =============================================================================
