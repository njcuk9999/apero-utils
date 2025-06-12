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
        self.basename = ''

    def load_rdb(self):
        abspath = os.path.join(self.path, self.basename)

        if not os.path.exists(abspath):
            raise FileNotFoundError(abspath)

        return Table.read(abspath)

# =============================================================================
# Define functions
# =============================================================================



# =============================================================================
# Define files
# =============================================================================
apero = File('apero')
apero.path = '/data/cook/nirps_comp/lbl/apero/lblrdb'
apero.basename = 'lbl_GL699_GL699.rdb'

eso = File('eso')
eso.path = '/data/cook/nirps_comp/lbl/eso/lblrdb'
eso.basename = 'lbl_GL699_GL699.rdb'


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # load tables
    apero_table = apero.load_rdb()
    eso_table = eso.load_rdb()

    apero_rjds = apero_table['rjd']
    apero_rvs = apero_table['vrad']
    apero_ervs = apero_table['svrad']

    eso_rjds = eso_table['rjd']
    eso_rvs = eso_table['vrad']
    eso_ervs = eso_table['svrad']

    # compute stats
    apero_sys_vel = np.nanmedian(apero_rvs)
    apero_p50 = np.nanpercentile(apero_ervs, 50)
    eso_sys_vel = np.nanmedian(eso_rvs)
    eso_p50 = np.nanpercentile(eso_ervs, 50)


    plt.errorbar(apero_rjds, apero_rvs, yerr=apero_ervs,
                 label=f'APERO ({apero_sys_vel:.2f} $\pm$ {apero_p50:.2f} m/s)',
                 marker='o', ls='None')
    plt.errorbar(eso_rjds, eso_rvs, yerr=eso_ervs,
                 label=f'ESO ({eso_sys_vel:.2f} $\pm$ {eso_p50:.2f} m/s)',
                 marker='o', ls='None')

    plt.legend(loc=0)

    plt.xlabel('MJD')
    plt.ylabel('LBL RV [m/s]')
    plt.savefig(os.path.join(PLOT_PATH, 'lbl_comp_raw.png'))
    plt.show()

    plt.errorbar(apero_rjds, apero_rvs - np.nanmedian(apero_rvs),
                 yerr=apero_ervs,
                 label=f'APERO ({apero_sys_vel:.2f} $\pm$ {apero_p50:.2f} m/s)',
                 marker='o', ls='None')
    plt.errorbar(eso_rjds, eso_rvs - np.nanmedian(eso_rvs),
                 yerr=eso_ervs,
                 label=f'ESO ({eso_sys_vel:.2f} $\pm$ {eso_p50:.2f} m/s)',
                 marker='o', ls='None')

    plt.legend(loc=0)

    plt.xlabel('MJD')
    plt.ylabel('LBL RV - Median(RV) [m/s]')
    plt.savefig(os.path.join(PLOT_PATH, 'lbl_comp_rm'
                                        '-med.png'))
    plt.show()



    # -------------------------------------------------------------------------



# =============================================================================
# End of code
# =============================================================================
