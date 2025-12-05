#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-06-10 at 12:31

@author: cook
"""
import os
import uu

import numpy as np
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from astropy.io import fits
from astropy.table import Table
from astropy import constants
from astropy import units as uu
from tqdm import tqdm
import warnings

# =============================================================================
# Define variables
# =============================================================================
BASENAME = 'NIRPS.2023-03-03T09:53:21.430'

ORDER_NUM = 60

speed_of_light_ms = constants.c.value

# Set global font size
mpl.rcParams.update({'font.size': 24})
# Set background color in plots
PLOT_BACKGROUND_COLOR = '#FEFDE1'
PLOT_BACKGROUND_COLOR2 = '#E1F0FE'

# set plot path
PLOT_PATH = '/data/cook/nirps_comp/plots'

MJD_KEY = 'MJD-OBS'

# =============================================================================
# Define classes
# =============================================================================
class File:
    def __init__(self, name):
        self.name = name
        self.path = ''
        self.obsdir = None
        self.basename = ''
        self.prefix = ''
        self.suffix = ''
        self.rv_key = None
        self.erv_key = None
        self.rv_units = None
        self.erv_units = None

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

    def load_header(self, **kwargs):
        abspath = self.get_abspath()
        return fits.getheader(abspath, **kwargs)

    def load_rvs(self, ext=0):
        # get all files
        files = self.get_files()
        # print progress
        print(f'Loading RVS for {len(files)} files [{self.name}].')
        # storage for return
        mjds, rvs, ervs = [], [], []
        # filename
        for filename in tqdm(files):
            hdr = fits.getheader(filename, ext=ext)

            rv = (float(hdr[self.rv_key]) * self.rv_units).to(uu.m/uu.s).value
            erv = (float(hdr[self.erv_key]) * self.erv_units).to(uu.m/uu.s).value
            mjds.append(float(hdr[MJD_KEY]))
            rvs.append(rv)
            ervs.append(erv)

        mjds = np.array(mjds)
        rvs = np.array(rvs)
        ervs = np.array(ervs)

        return mjds, rvs, ervs


# =============================================================================
# Define functions
# =============================================================================



# =============================================================================
# Define files
# =============================================================================
apero = File('apero')
apero.path = '/cosmos99/nirps/apero-data/nirps_he_online/objects/GL699'
apero.obsdir = None
apero.basename = '*'
apero.suffix = 'v.fits'
apero.rv_key = 'RV_OBJ'
apero.erv_key = 'DVRMS_CC'
apero.rv_units = uu.km/uu.s
apero.erv_units = uu.m/uu.s

eso = File('eso')
eso.path = '/data/cook/nirps_comp/geneva/GL699/'
eso.obsdir = '*'
eso.basename = '*'
eso.prefix = 'r.'
eso.suffix = '_CCF_TELL_CORR_A.fits'
eso.rv_key = 'HIERARCH ESO QC CCF RV'
eso.erv_key = 'HIERARCH ESO QC CCF RV ERROR'
eso.rv_units = uu.km/uu.s
eso.erv_units = uu.km/uu.s

# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # load rvs
    apero_mjds, apero_rvs, apero_ervs = apero.load_rvs(ext=1)
    eso_mjds, eso_rvs, eso_ervs = eso.load_rvs()

    # compute stats
    apero_sys_vel = np.nanmedian(apero_rvs)
    apero_p50 = np.nanpercentile(apero_ervs, 50)
    eso_sys_vel = np.nanmedian(eso_rvs)
    eso_p50 = np.nanpercentile(eso_ervs, 50)

    plt.errorbar(apero_mjds, apero_rvs, yerr=apero_ervs,
                 label=f'APERO ({apero_sys_vel:.2f} $\pm$ {apero_p50:.2f} m/s)',
                 marker='o', ls='None')
    plt.errorbar(eso_mjds, eso_rvs, yerr=eso_ervs,
                 label=f'ESO ({eso_sys_vel:.2f} $\pm$ {eso_p50:.2f} m/s)',
                 marker='o', ls='None')

    plt.legend(loc=0)

    plt.xlabel('MJD')
    plt.ylabel('CCF RV [m/s]')
    plt.savefig(os.path.join(PLOT_PATH, 'ccf_comp_raw.png'))
    plt.show()

    plt.errorbar(apero_mjds, apero_rvs - np.nanmedian(apero_rvs),
                 yerr=apero_ervs,
                 label=f'APERO ({apero_sys_vel:.2f} $\pm$ {apero_p50:.2f} m/s)',
                 marker='o', ls='None')
    plt.errorbar(eso_mjds, eso_rvs - np.nanmedian(eso_rvs),
                 yerr=eso_ervs,
                 label=f'ESO ({eso_sys_vel:.2f} $\pm$ {eso_p50:.2f} m/s)',
                 marker='o', ls='None')

    plt.legend(loc=0)

    plt.xlabel('MJD')
    plt.ylabel('CCF RV - Median(RV) [m/s]')
    plt.savefig(os.path.join(PLOT_PATH, 'ccf_comp_rm'
                                        '-med.png'))
    plt.show()



    # -------------------------------------------------------------------------



# =============================================================================
# End of code
# =============================================================================
