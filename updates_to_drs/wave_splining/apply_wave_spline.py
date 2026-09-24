#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
apply_wave_spline.py

Find all *wave* FITS files in WAVEDIR and up-sample each one along the
pixel axis from WAVESIZE1 pixels to WAVESIZE2 pixels using a simple
half-pixel gradient interpolation, overwriting the file in place.

Created on 2026-09-22

@author: cook
"""
import glob
import os

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

# =============================================================================
# Define variables
# =============================================================================
# Directory to search for wave files
WAVEDIR = '/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero-drs/apero/apero-assets/spirou/calib/'
# WAVEDIR = '/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero-drs/apero/apero-assets/nirps_he/calib/'
# WAVEDIR = '/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero-drs/apero/apero-assets/nirps_ha/calib/'


# Glob pattern used to find wave files inside WAVEDIR
WAVEPAT = '*wave*.fits'

# Expected pixel size of the input wave array (columns)
WAVESIZE1 = 4088

# Target pixel size of the output wave array (columns)
WAVESIZE2 = 8176

MODE = 'MAKE'
# MODE = 'TEST'


# =============================================================================
# Define functions
# =============================================================================
def find_wave_files(directory: str, pattern: str) -> list:
    """
    Return a sorted list of wave file paths matching pattern in directory.

    :param directory: str, directory to search
    :param pattern: str, glob pattern (e.g. '*wave*.fits')
    :return: list of str, absolute file paths found
    """
    search_path = os.path.join(directory, pattern)
    matches = sorted(glob.glob(search_path))
    return matches


def find_wave_extension(hdul: fits.HDUList) -> int:
    """
    Return the index of the HDU that contains the wavelength array.

    Checks for a named 'WAVESOL_REF' ImageHDU first (NIRPS layout); falls
    back to extension 0 (SPIROU / generic PrimaryHDU layout).

    :param hdul: fits.HDUList, open FITS file
    :return: int, HDU index that holds the wavelength data
    """
    for idx, hdu in enumerate(hdul):
        if hdu.name == 'WAVESOL_REF':
            return idx
    # default to primary extension
    return 0


def make_wave_spline(wavefile: str) -> None:
    """
    Up-sample a wave file from WAVESIZE1 to WAVESIZE2 pixels, overwriting
    the file in place.  Every even column keeps the original value; every
    odd column is set to original + gradient/2 so the spacing is uniform
    at half the original pixel scale.

    Supports both single-extension files (SPIROU: data in PRIMARY) and
    multi-extension files (NIRPS: wavelength data in 'WAVESOL_REF' ImageHDU,
    with additional extensions preserved unchanged).

    :param wavefile: str, path to the FITS wave file to process
    :return: None
    """
    with fits.open(wavefile) as hdul:
        # locate the extension that holds the wavelength array
        ext = find_wave_extension(hdul)
        wave1 = hdul[ext].data
        # validate shape matches expected input size
        shape1 = wave1.shape
        if shape1[1] != WAVESIZE1:
            emsg = ('Wave file {0} has pixel axis size {1} '
                    'but expected {2} — skipping.')
            print(emsg.format(wavefile, shape1[1], WAVESIZE1))
            return
        # deep-copy all HDUs so we can modify one and write the full list
        hdul_out = fits.HDUList(
            [hdu.copy() for hdu in hdul]
        )

    # allocate output array with double the pixel columns
    wave2 = np.zeros([shape1[0], WAVESIZE2])
    # even columns: copy original values
    wave2[:, 0::2] = wave1
    # odd columns: original + half-pixel gradient offset
    wave2[:, 1::2] = wave1 + np.gradient(wave1, axis=1) / 2.0

    # replace only the wavelength extension data; all other HDUs are unchanged
    hdul_out[ext].data = wave2

    print(f'Writing up-sampled wave to {wavefile} (ext={ext})')
    hdul_out.writeto(wavefile, overwrite=True)


def test_wave_spline(wavefile: str) -> None:
    """
    Visual test: overlay the original (reconstructed from every-other
    column) and up-sampled wave grids for a single file.

    Works with both single-extension (SPIROU) and multi-extension (NIRPS)
    files by delegating extension detection to find_wave_extension.

    :param wavefile: str, path to the already-processed FITS wave file
    :return: None
    """
    with fits.open(wavefile) as hdul:
        ext = find_wave_extension(hdul)
        wave2 = hdul[ext].data

    # reconstruct the original-resolution grid from even columns
    wave1 = wave2[:, 0::2]
    # build matching fake spectra
    data1 = np.arange(0, np.prod(wave1.shape), 1).reshape(wave1.shape)
    data2 = np.arange(0, np.prod(wave1.shape), 0.5).reshape(wave2.shape)

    plot_spec(wave1, data1, label='Original',
              colours=['orange', 'purple'], marker='x')
    plot_spec(wave2, data2, label='Up-sampled',
              colours=['red', 'green'], marker='+')
    plt.title(os.path.basename(wavefile))
    plt.legend()
    plt.show()


def plot_spec(wave: np.ndarray, spectrum: np.ndarray,
              label: str, colours: list, marker: str) -> None:
    """
    Plot every order of a 2-D spectrum against its wave grid.

    :param wave: np.ndarray, shape (n_orders, n_pixels), wavelength grid
    :param spectrum: np.ndarray, shape (n_orders, n_pixels), spectral data
    :param label: str, legend label (applied only to the first two orders)
    :param colours: list of two str, alternating colours for even/odd orders
    :param marker: str, matplotlib marker style
    :return: None
    """
    for order_num in range(spectrum.shape[0]):
        colour = colours[0] if order_num % 2 == 0 else colours[1]
        # only attach the label to the first two orders to avoid legend clutter
        order_label = label if order_num in [0, 1] else None
        plt.plot(wave[order_num], spectrum[order_num],
                 color=colour, label=order_label, marker=marker)


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":

    # locate all wave files in the target directory
    wave_files = find_wave_files(WAVEDIR, WAVEPAT)

    if not wave_files:
        print(f'No files matching "{WAVEPAT}" found in {WAVEDIR}')
    else:
        print(f'Found {len(wave_files)} wave file(s) in {WAVEDIR}')

    if MODE == 'MAKE':
        for wfile in wave_files:
            make_wave_spline(wfile)
    elif MODE == 'TEST':
        for wfile in wave_files:
            test_wave_spline(wfile)
    else:
        raise ValueError(f'Unknown MODE: {MODE}')


# =============================================================================
# End of code
# =============================================================================
