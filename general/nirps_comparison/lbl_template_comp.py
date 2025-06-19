#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-06-17 at 14:27

@author: cook
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from tqdm import tqdm
import warnings
import os
from scipy.interpolate import InterpolatedUnivariateSpline

# =============================================================================
# Define variables
# =============================================================================

PATH = '/data/cook/nirps_comp/lbl'
PATH = '/scratch2/nirps/misc/nirps_comp/lbl/'

apero_template = os.path.join(PATH, 'apero/templates/Template_s1dv_GL699_sc1d_v_file_A.fits')
eso_template = os.path.join(PATH, 'eso/templates/Template_GL699_NIRPS_HE_ESO.fits')

apero_mask = os.path.join(PATH, 'apero/masks/GL699_full.fits')
eso_mask = os.path.join(PATH, 'eso/masks/GL699_full.fits')
# -----------------------------------------------------------------------------
LBL_KEY = 'SYSTVELO'

# FLUX_COL_APERO = 'flux_odd'
# FLUX_COL_ESO = 'flux_even'

FLUX_COL_APERO = 'flux'
FLUX_COL_ESO = 'flux'

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
    # load the data
    apero_table = Table.read(apero_template)
    eso_table = Table.read(eso_template)
    apero_hdr = fits.getheader(apero_mask)
    eso_hdr = fits.getheader(eso_mask)

    apero_wave = np.array(apero_table['wavelength']).astype(float)
    apero_flux = np.array(apero_table[FLUX_COL_APERO]).astype(float)
    eso_wave = np.array(eso_table['wavelength']).astype(float)
    eso_flux = np.array(eso_table[FLUX_COL_ESO]).astype(float)



    eso_mask = np.isfinite(eso_wave) & np.isfinite(eso_flux)

    eso_spline = InterpolatedUnivariateSpline(eso_wave[eso_mask], eso_flux[eso_mask])
    diff = apero_flux - eso_spline(apero_wave)
    d16, d50, d84 = np.nanpercentile(diff, [16, 50, 84])
    dsig = (d84 - d16) / 2

    ratio = apero_flux / eso_spline(apero_wave)
    r16, r50, r84 = np.nanpercentile(ratio, [16, 50, 84])
    rsig = (r84 - r16) / 2

    # get sys velo
    apero_vsys = float(apero_hdr[LBL_KEY])
    eso_vsys = float(eso_hdr[LBL_KEY])
    # plot the apero data
    plt.close()
    fig, frames = plt.subplots(nrows=3, ncols=1, figsize=(10, 6),
                               sharex='all')

    frames[0].plot(apero_wave, apero_flux, label=f'APERO SYSVELO={apero_vsys:.3f} m/s')
    frames[0].plot(eso_wave, eso_flux, label=f'ESO SYSTVELO={eso_vsys:.3f} m/s')

    frames[1].plot(apero_wave, diff, color='k')
    frames[1].set(ylim=[d50 - 5 * dsig, d50 + 5* dsig],
                  xlabel='Wavelength [nm]', ylabel='APERO - ESO')

    frames[2].plot(apero_wave, ratio, color='k')
    frames[2].set(ylim=[r50 - 5 * rsig, r50 + 5* rsig],
                  xlabel='Wavelength [nm]', ylabel='APERO / ESO')

    frames[0].legend(loc=0)

    plt.suptitle('GL699 LBL Template')


    plt.show()

# =============================================================================
# End of code
# =============================================================================
