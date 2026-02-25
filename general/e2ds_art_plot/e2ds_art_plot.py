#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-03 at 14:51

@author: cook
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# =============================================================================
# Define variables
# =============================================================================
IPARAMS = dict()
IPARAMS['ilocater'] = dict()
IPARAMS['ilocater']['WAVE'] = '/scratch2/ilocater/drs-data/ilocater_test_2/calib/cooldown13_composite_fp_hxrg_pp_e2dsff_A_wave_night_A.fits'
IPARAMS['ilocater']['E2DS'] = '/scratch2/ilocater/drs-data/ilocater_test_2/red/2026-01-30/iLocater_lab_20260130_0055_hxrgproc_pp_e2dsff_A.fits'
IPARAMS['ilocater']['NAME'] = 4

IPARAMS['spirou'] = dict()
IPARAMS['spirou']['WAVE'] = '/scratch2/spirou/drs-data/spirou_xxs_07/calib/0F66E13BEDa_pp_e2dsff_AB_wave_night_AB.fits'
IPARAMS['spirou']['E2DS'] = '/scratch2/spirou/drs-data/spirou_xxs_07/red/2020-08-31/2510301o_pp_e2dsff_AB.fits'
IPARAMS['spirou']['NAME'] = 1

IPARAMS['nirps-he'] = dict()
IPARAMS['nirps-he']['WAVE'] = '/scratch2/nirps/misc/nirps_proxima/C1CAAEE1D9_pp_e2dsff_A_wave_night_A.fits'
IPARAMS['nirps-he']['E2DS'] = '/scratch2/nirps/misc/nirps_proxima/NIRPS_2024-09-04T00_09_45_760_pp_e2dsff_A.fits'
IPARAMS['nirps-he']['NAME'] = 3

IPARAMS['nirps_ha'] = dict()
IPARAMS['nirps_ha']['WAVE'] = ''
IPARAMS['nirps_ha']['E2DS'] = ''
IPARAMS['nirps_ha']['NAME'] = None

# -----------------------------------------------------------------------------
OUTDIR = '/scratch2/ilocater/misc/art_plot/'



# =============================================================================
# Define functions
# =============================================================================
def plot_e2ds(instrument):

    # Get parameters for the instrument
    wave_file = IPARAMS[instrument]['WAVE']
    e2ds_file = IPARAMS[instrument]['E2DS']
    name = IPARAMS[instrument]['NAME']

    if not os.path.exists(wave_file):
        print(f"Wave file for {instrument} not found: {wave_file}")
        return
    if not os.path.exists(e2ds_file):
        print(f"E2DS file for {instrument} not found: {e2ds_file}")
        return

    flux = fits.getdata(e2ds_file)
    wave = fits.getdata(wave_file)
    # Scale factor for y-axis (number of pixels per order)
    ypix_scale = int(flux.shape[1] / flux.shape[0] / 2)
    gap_scale = int(ypix_scale / 10)  # Number of pixels gap between orders

    fig, frame = plt.subplots(ncols=1, nrows=1, figsize=(20, 20))
    # turn off axis
    frame.set_axis_off()
    # set background colour to black
    frame.set_facecolor('black')
    fig.patch.set_facecolor('black')

    # build an image that is scaled in the y direction by YPIX_SCALE
    # with gaps between orders
    total_height = flux.shape[0] * (ypix_scale + gap_scale)
    flux_image = np.zeros((total_height, flux.shape[1]))
    wave_image = np.zeros((total_height, wave.shape[1]))

    for i in range(flux.shape[0]):
        start_row = i * (ypix_scale + gap_scale)
        end_row = start_row + ypix_scale
        flux_image[start_row:end_row, :] = flux[i, :]
        wave_image[start_row:end_row, :] = wave[i, :]

    # Get min/max wavelength from the actual data (not including zeros)
    wave_min = np.nanmin(wave)
    wave_max = np.nanmax(wave)

    # Normalize flux to 0-1 range for alpha (transparency)
    flux_min, flux_max = np.nanmin(flux), np.nanmax(flux)
    flux_norm = (flux_image - flux_min) / (flux_max - flux_min)
    # Handle any NaN or invalid values
    flux_norm = np.nan_to_num(flux_norm, nan=0.0)

    # Use wavelength for color and flux for transparency
    # Set vmin/vmax to use full colormap range
    plt.imshow(wave_image, aspect='auto', origin='lower', cmap='jet',
               alpha=flux_norm, vmin=wave_min, vmax=wave_max,
               interpolation='None')

    plt.savefig(OUTDIR + f'e2ds_art_{name}.png', dpi=300, bbox_inches='tight')
    plt.close()

# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------:
    # run main
    for _instrument in IPARAMS:
        print(f'Processing instrument: {_instrument}')
        plot_e2ds(_instrument)

# =============================================================================
# End of code
# =============================================================================

