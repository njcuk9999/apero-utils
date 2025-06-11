#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-06-10 at 14:41

@author: cook
"""
import os
import glob
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
APERO_BERV_KEY = 'BERV'
ESO_BERV_KEY = 'HIERARCH ESO QC BERV'
MJD_KEY = 'MJD-OBS'

# =============================================================================
# Define functions
# =============================================================================
class File:
    def __init__(self):
        self.path = ''
        self.obsdir = None
        self.basename = ''
        self.prefix = ''
        self.suffix = ''
        self.format = None

    def get_abspath(self, check: bool = True):
        basename = self.prefix + self.basename + self.suffix
        if self.obsdir is None:
            abspath = os.path.join(self.path,basename)
        else:
            abspath = os.path.join(self.path, self.obsdir, basename)
        if not os.path.exists(abspath) and check:
            raise ValueError(f'File {abspath} does not exist.')
        return abspath

    def get_files(self):
        return glob.glob(self.get_abspath(check=False))

    def load_data(self, **kwargs):
        """Load the file data."""
        abspath = self.get_abspath()

        if self.format == 'image':
            return fits.getdata(abspath, **kwargs)
        elif self.format == 'table':
            return Table.read(abspath, **kwargs)

    def load_header(self, **kwargs):
        abspath = self.get_abspath()
        return fits.getheader(abspath, **kwargs)

    def get_base(self, filename):

        basename = os.path.basename(filename)

        if len(self.prefix) > 0:
            basename = basename.split(self.prefix)[-1]

        if len(self.suffix) > 0:
            basename = basename.split(self.suffix)[0]

        return basename

# =============================================================================
# Define files
# =============================================================================

apero = File()
apero.path = '/cosmos99/nirps/apero-data/nirps_he_online/objects/GL699'
apero.obsdir = None
apero.basename = '*'
apero.suffix = 'e.fits'
apero.format = 'image'

eso = File()
eso.path = '/data/cook/nirps_comp/geneva/GL699/'
eso.obsdir = '*'
eso.basename = '*'
eso.prefix = 'r.'
eso.suffix = '_S2D_BLAZE_A.fits'
eso.format = 'image'



# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # get list of apero files
    apero_files = apero.get_files()
    # get list of eso files
    eso_files = eso.get_files()

    apero_bervs = dict()
    apero_mjds = dict()
    eso_bervs = dict()
    eso_mjds = dict()


    print('Loading BERVs from APERO files...')
    for apero_filename in tqdm(apero_files):

        basename = apero.get_base(apero_filename)
        hdr = fits.getheader(apero_filename, ext=1)
        berv = float(hdr[APERO_BERV_KEY])
        mjd = hdr[MJD_KEY]

        apero_bervs[basename] = berv
        apero_mjds[basename] = mjd

    print('Loading BERVs from ESO files...')
    for eso_filename in tqdm(eso_files):
        basename = eso.get_base(eso_filename)
        hdr = fits.getheader(eso_filename)
        berv = float(hdr[ESO_BERV_KEY])
        mjd = hdr[MJD_KEY]

        eso_bervs[basename] = berv
        eso_mjds[basename] = mjd


    # -------------------------------------------------------------------------


    # plot points that are in common between both eso and apero
    common_bases = set(apero_bervs.keys()).intersection(set(eso_bervs.keys()))

    # make a list of values for common bases
    c_apero_bervs, c_apero_mjds = [], []
    for base in common_bases:
        c_apero_bervs.append(apero_bervs[base])
        c_apero_mjds.append(apero_mjds[base])
    c_eso_bervs, c_eso_mjds = [], []
    for base in common_bases:
        c_eso_bervs.append(eso_bervs[base])
        c_eso_mjds.append(eso_mjds[base])
    # make a list of values not in other set
    u_apero_bervs, u_apero_mjds = [], []
    for base in apero_bervs.keys():
        if base not in common_bases:
            u_apero_bervs.append(apero_bervs[base])
            u_apero_mjds.append(apero_mjds[base])
    u_eso_bervs, u_eso_mjds = [], []
    for base in eso_bervs.keys():
        if base not in common_bases:
            u_eso_bervs.append(eso_bervs[base])
            u_eso_mjds.append(eso_mjds[base])

    # push to numpy arrays
    c_apero_bervs = np.array(c_apero_bervs)
    c_apero_mjds = np.array(c_apero_mjds)
    c_eso_bervs = np.array(c_eso_bervs)
    c_eso_mjds = np.array(c_eso_mjds)
    u_apero_bervs = np.array(u_apero_bervs)
    u_apero_mjds = np.array(u_apero_mjds)
    u_eso_bervs = np.array(u_eso_bervs)
    u_eso_mjds = np.array(u_eso_mjds)


    # plot berv
    plt.close()
    fig, frames = plt.subplots(figsize=(10, 6), nrows=2, ncols=1, sharex='all')

    frames[0].plot(c_apero_mjds, 1000 * c_apero_bervs, 'o', color='orange',
                   label='APERO BERVs')
    frames[0].plot(c_eso_mjds, 1000 * c_eso_bervs, 'o', color='blue',
                   label='ESO BERVs')

    frames[0].plot(u_apero_mjds, 1000 * u_apero_bervs, 'x', color='orange',
                   label=f'APERO BERVs (unmatched N={len(u_apero_mjds)})')
    frames[0].plot(u_eso_mjds, 1000 * u_eso_bervs, 'x', color='blue',
                   label=f'ESO BERVs (unmatched N={len(u_eso_mjds)})')
    frames[0].set(xlabel='MJD', ylabel='BERV (m/s)')

    # frame 1 is residuals
    frames[1].plot(c_apero_mjds, 1000 * (c_apero_bervs - c_eso_bervs),
                   'o', color='green', label='APERO - ESO')
    frames[1].set(xlabel='MJD', ylabel='BERV residuals (m/s)')


    frames[0].legend(loc=0)
    frames[1].legend(loc=0)
    frames[0].grid(True, which='major', linestyle='-', alpha=0.5)
    frames[0].grid(True, which='minor', linestyle='--', alpha=0.25)
    frames[1].grid(True, which='major', linestyle='-', alpha=0.5)
    frames[1].grid(True, which='minor', linestyle='--', alpha=0.25)

    plt.suptitle('GL699 BERV Comparison')

    plt.tight_layout()

    plt.show()

    # sort them
    apero_sort_mask = np.argsort(c_apero_mjds)
    eso_sort_mask = np.argsort(c_eso_mjds)

    c_apero_mjds = c_apero_mjds[apero_sort_mask]
    c_apero_bervs = c_apero_bervs[apero_sort_mask]
    c_eso_mjds = c_eso_mjds[eso_sort_mask]
    c_eso_bervs = c_eso_bervs[eso_sort_mask]

    mjd_min = 60460
    mjd_max = 60461

    apero_mask = (c_apero_mjds > mjd_min) & (c_apero_mjds < mjd_max)
    eso_mask = (c_eso_mjds > mjd_min) & (c_eso_mjds < mjd_max)

    apero_coeffs = np.polyfit(c_apero_mjds[apero_mask],
                              1000 * c_apero_bervs[apero_mask],
                              1)
    eso_coeffs = np.polyfit(c_eso_mjds[eso_mask],
                            1000 * c_eso_bervs[eso_mask],
                            1)
    mjd_grid = np.arange(mjd_min, mjd_max, 0.01)

    plt.close()
    fig, frames = plt.subplots(figsize=(10, 6), nrows=2, ncols=1, sharex='all')

    frames[0].plot(c_apero_mjds[apero_mask], 1000 * c_apero_bervs[apero_mask],
             marker='o', ls='None', color='orange', label='APERO BERVs')

    frames[0].plot(c_apero_mjds[apero_mask],
             np.polyval(apero_coeffs, c_apero_mjds[apero_mask]),
             color='orange', label='APERO linear fit')

    frames[0].plot(c_eso_mjds[eso_mask], 1000 * c_eso_bervs[eso_mask],
             marker='o', ls='None', color='blue', label='ESO BERVs')

    frames[0].plot(c_eso_mjds[eso_mask],
             np.polyval(eso_coeffs, c_eso_mjds[eso_mask]),
             color='blue', label='ESO linear fit')

    frames[1].plot(c_apero_mjds[apero_mask], 1000 * c_apero_bervs[apero_mask] - np.polyval(apero_coeffs, c_apero_mjds[apero_mask]),
             marker='o', ls='None', color='orange', label='APERO BERVs - APERO linear fit')
    frames[1].plot(c_eso_mjds[eso_mask], 1000 * c_eso_bervs[eso_mask] - np.polyval(eso_coeffs, c_eso_mjds[eso_mask]),
             marker='o', ls='None', color='blue', label='ESO BERVs - ESO linear fit')

    plt.legend(loc=0)

    plt.show()

# =============================================================================
# End of code
# =============================================================================
