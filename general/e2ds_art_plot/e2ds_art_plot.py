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
from astropy.table import Table
from astropy.visualization import ImageNormalize, LinearStretch, ZScaleInterval

# =============================================================================
# Define variables
# =============================================================================
IPARAMS = dict()
IPARAMS['ilocater'] = dict()
IPARAMS['ilocater']['WAVE'] = '/scratch2/ilocater/drs-data/ilocater_test_2/calib/cooldown13_composite_fp_hxrg_pp_e2dsff_A_wave_night_A.fits'
IPARAMS['ilocater']['E2DS'] = '/scratch2/ilocater/drs-data/ilocater_test_2/red/2026-01-30/iLocater_lab_20260130_0055_hxrgproc_pp_e2dsff_A.fits'
IPARAMS['ilocater']['NAME'] = 'iLocater'
IPARAMS['ilocater']['S1D'] = '/scratch2/ilocater/drs-data/ilocater_test_2/red/2026-01-30/iLocater_lab_20260130_0055_hxrgproc_pp_s1d_v_A.fits'
IPARAMS['ilocater']['Title'] = 'iLocater [Extracted] In Lab: Sun'

IPARAMS['spirou'] = dict()
IPARAMS['spirou']['WAVE'] = '/scratch2/spirou/drs-data/spirou_xxs_07/calib/0F66E13BEDa_pp_e2dsff_AB_wave_night_AB.fits'
IPARAMS['spirou']['E2DS'] = '/scratch2/spirou/drs-data/spirou_xxs_07/red/2020-08-31/2510301o_pp_e2dsff_AB.fits'
IPARAMS['spirou']['NAME'] = 'SPIRou'
IPARAMS['spirou']['s.fits'] = '/scratch2/spirou/drs-data/spirou_xxs_07/out/2020-08-31/2510301s.fits'
IPARAMS['spirou']['Title'] = 'SPIRou [Telluric Corrected] GL699'

IPARAMS['nirps-he'] = dict()
IPARAMS['nirps-he']['WAVE'] = '/scratch2/nirps/misc/nirps_proxima/C1CAAEE1D9_pp_e2dsff_A_wave_night_A.fits'
IPARAMS['nirps-he']['E2DS'] = '/scratch2/nirps/misc/nirps_proxima/NIRPS_2024-09-04T00_09_45_760_pp_e2dsff_A.fits'
IPARAMS['nirps-he']['NAME'] = 'NIRPS-HE'
IPARAMS['nirps-he']['s.fits'] = '/scratch2/nirps/misc/nirps_proxima/NIRPS.2024-09-04T00:09:45.760s.fits'
IPARAMS['nirps-he']['Title'] = 'NIRPS-HE [Telluric Corrected] PROXIMA'

IPARAMS['spip'] = dict()
IPARAMS['spip']['WAVE'] = '/scratch2/spip/drs-data/spip_mini_07/calib/8D18FA61BCa_pp_e2dsff_AB_wave_night_AB.fits'
IPARAMS['spip']['E2DS'] = '/scratch2/spip/drs-data/spip_mini_07/red/2026-03-02/1F2749A0BAc_pp_e2dsff_AB.fits'
IPARAMS['spip']['NAME'] = 'SPIP'
IPARAMS['spip']['S1D'] = '/scratch2/spip/drs-data/spip_mini_07/red/2026-03-02/1F2749A0BAc_pp_s1d_v_AB.fits'
IPARAMS['spip']['Title'] = 'SPIP [Extracted] UNe Hollow-Cathode'

# -----------------------------------------------------------------------------
OUTDIR = '/scratch2/apero/casca_2026_poster/e2ds_art/'

PLOT_S1D = True
# S1D mid point 0 = bottom 1 = top
S1D_Y_CENT = 0.20
# S1D extent as a fraction of the image height
S1D_Y_EXTENT = 0.30

CMAP = 'inferno'
# Minimum per-column finite fraction used to keep a column in wave stretch
WAVE_VALID_COLUMN_FRAC = 0.15
TITLE_FONT_SIZE = 28




# =============================================================================
# Define functions
# =============================================================================
def plot_s1d(instrument, wave_min, wave_max, frame, wave_image):
    if 'S1D' in IPARAMS[instrument]:
        s1d_table = Table.read(IPARAMS[instrument]['S1D'],
                               hdu='EXT_S1D_V')
        s1d_wave, s1d_flux = s1d_table['wavelength'], s1d_table['flux']
    elif 's.fits' in IPARAMS[instrument]:
        s1d_table = Table.read(IPARAMS[instrument]['s.fits'],
                               hdu='UniformVelocity')
        s1d_wave, s1d_flux = s1d_table['Wave'], s1d_table['FluxATelluCorrected']

    else:
        return frame

    s1d_wave = np.asarray(s1d_wave, dtype=float)
    s1d_flux = np.asarray(s1d_flux, dtype=float)

    # 1) Trim only contiguous invalid regions at the edges.
    joint_valid = np.isfinite(s1d_wave) & np.isfinite(s1d_flux)
    if not np.any(joint_valid):
        return frame
    valid_idx = np.flatnonzero(joint_valid)
    s1d_wave = s1d_wave[valid_idx[0]:valid_idx[-1] + 1]
    s1d_flux = s1d_flux[valid_idx[0]:valid_idx[-1] + 1]

    # 3) Median + 5 sigma outlier mask in flux.
    finite_flux = s1d_flux[np.isfinite(s1d_flux)]
    if finite_flux.size == 0:
        return frame
    flux_med = np.nanmedian(finite_flux)
    flux_sig = np.nanstd(finite_flux)
    if np.isfinite(flux_sig) and flux_sig > 0:
        high_mask = s1d_flux > (flux_med + 5.0 * flux_sig)
        s1d_flux = s1d_flux.copy()
        s1d_flux[high_mask] = np.nan

    # 2) Scale wavelength/flux to fit the image overlay footprint.
    finite_wave = s1d_wave[np.isfinite(s1d_wave)]
    if finite_wave.size == 0:
        return frame
    s1d_wave_min = np.nanmin(finite_wave)
    s1d_wave_max = np.nanmax(finite_wave)
    wave_span = s1d_wave_max - s1d_wave_min
    if np.isfinite(wave_span) and wave_span > 0:
        s1d_wave_norm = (s1d_wave - s1d_wave_min) / wave_span
    else:
        # Fallback to image wave stretch if S1D wavelength is degenerate.
        fallback_span = wave_max - wave_min
        if not np.isfinite(fallback_span) or fallback_span <= 0:
            return frame
        s1d_wave_norm = (s1d_wave - wave_min) / fallback_span

    s1d_wave_norm = np.asarray(s1d_wave_norm, dtype=float)
    s1d_wave_norm[(s1d_wave_norm < 0.0) | (s1d_wave_norm > 1.0)] = np.nan

    finite_flux = s1d_flux[np.isfinite(s1d_flux)]
    if finite_flux.size == 0:
        return frame
    flux_min = np.nanmin(finite_flux)
    flux_span = np.nanmax(finite_flux) - flux_min
    if np.isfinite(flux_span) and flux_span > 0:
        s1d_flux_norm = (s1d_flux - flux_min) / flux_span
    else:
        s1d_flux_norm = np.full_like(s1d_flux, 0.5, dtype=float)
    s1d_flux_norm = np.asarray(s1d_flux_norm, dtype=float)
    s1d_flux_norm[~np.isfinite(s1d_flux)] = np.nan

    image_width = wave_image.shape[1]
    image_height = wave_image.shape[0]
    s1d_x = s1d_wave_norm * (image_width - 1)
    s1d_y_center = S1D_Y_CENT * (image_height - 1)
    s1d_y_span = S1D_Y_EXTENT * (image_height - 1)
    s1d_y = s1d_y_center + (s1d_flux_norm - 0.5) * s1d_y_span



    frame.plot(s1d_x, s1d_y, color='white', linewidth=1.0, zorder=10,
               alpha=0.75)

    return frame


def get_wave_stretch_limits(wave, min_valid_frac=WAVE_VALID_COLUMN_FRAC):
    # Collapse 2D wave map by column and trim NaN-heavy edge regions
    finite_mask = np.isfinite(wave)
    valid_frac_by_col = np.mean(finite_mask, axis=0)
    good_cols = valid_frac_by_col >= min_valid_frac

    if np.any(good_cols):
        good_idx = np.flatnonzero(good_cols)
        split_points = np.where(np.diff(good_idx) > 1)[0] + 1
        runs = np.split(good_idx, split_points)
        largest_run = max(runs, key=len)
        first_col = int(largest_run[0])
        last_col = int(largest_run[-1])
        wave_trim = wave[:, first_col:last_col + 1]
        if np.isfinite(wave_trim).any():
            return np.nanmin(wave_trim), np.nanmax(wave_trim)

    # Fallback to all finite data if no columns pass the threshold
    if np.isfinite(wave).any():
        return np.nanmin(wave), np.nanmax(wave)

    return np.nan, np.nan

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

    # Get min/max wavelength after trimming NaN-heavy edge columns
    wave_min, wave_max = get_wave_stretch_limits(wave)
    if not (np.isfinite(wave_min) and np.isfinite(wave_max)):
        print(f"Wave data for {instrument} has no finite values")
        plt.close()
        return

    # Normalize flux to 0-1 alpha values with a robust zscale + linear stretch
    finite_flux = flux[np.isfinite(flux)]
    if finite_flux.size > 0:
        zscale = ZScaleInterval()
        z_vmin, z_vmax = zscale.get_limits(finite_flux)
        flux_norm = ImageNormalize(vmin=z_vmin, vmax=z_vmax,
                                   stretch=LinearStretch(),
                                   clip=True)(flux_image)
    else:
        flux_norm = np.zeros_like(flux_image, dtype=float)
    # Keep gaps/invalid pixels fully transparent
    flux_norm = np.nan_to_num(flux_norm, nan=0.0)

    # Use wavelength for color and flux for transparency
    # Set vmin/vmax to use full colormap range
    frame.imshow(wave_image, aspect='auto', origin='lower', cmap=CMAP,
                 alpha=flux_norm, vmin=wave_min, vmax=wave_max,
                 interpolation='None')

    # plot the s1d
    if PLOT_S1D:
        plot_s1d(instrument, wave_min, wave_max, frame, wave_image)

    title_text = IPARAMS[instrument].get('Title', f"{name} E2DS + S1D")
    frame.text(0.5, 0.03, title_text, transform=frame.transAxes,
               color='white', fontsize=TITLE_FONT_SIZE, fontweight='bold',
               ha='center', va='top', wrap=False)

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

