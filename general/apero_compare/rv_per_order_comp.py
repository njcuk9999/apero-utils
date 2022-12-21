#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-12-21 at 15:59

@author: cook
"""
from astropy.table import Table
from astropy import constants
import os
import matplotlib.pyplot as plt


# =============================================================================
# Define variables
# =============================================================================
name1 = 'spirou@rali'
name2 = 'cook@jupiter'
name3 = 'cook@nb19'

ref_name = 'spirou@rali'

paths = dict()
paths[name1] = '/scratch2/spirou/drs-data/minidata2_2022-12-15_rali/red/2020-08-02'
paths[name2] = '/scratch2/spirou/drs-data/minidata2_07XXX/reduced/2020-08-02'
paths[name3] = '/scratch2/spirou/drs-data/minidata2_2022-12-15_nb19/red/2020-08-02'


filename = '2503782o_pp_e2dsff_tcorr_AB_ccf_gl699_neg.fits_AB.fits'

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
    data_ref = Table.read(os.path.join(paths[ref_name], filename), hdu=2)

    for name in paths:

        if name == ref_name:
            continue

        data_comp = Table.read(os.path.join(paths[name], filename), hdu=2)

        diff = 1000 * (data_comp['RV'] - data_ref['RV'])

        plt.plot(data_comp['Orders'], diff, label=f'{ref_name}-{name}',
                 marker='o', ls='None')

    plt.legend(loc=0)
    plt.xlabel('Orders')
    plt.ylabel('$\Delta$RV  m/s')
    plt.show()
    plt.close()

# =============================================================================
# End of code
# =============================================================================
