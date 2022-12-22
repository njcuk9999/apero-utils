#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-12-21 at 15:59

@author: cook
"""
from astropy.table import Table
import glob
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
paths[name1] = '/scratch2/spirou/drs-data/minidata2_2022-12-15_rali/red/'
paths[name2] = '/scratch2/spirou/drs-data/minidata2_07XXX/reduced/'
paths[name3] = '/scratch2/spirou/drs-data/minidata2_2022-12-15_nb19/red/'

obs_dir = '2020-10-07'

filestr = '*ccf_gl699_neg*_AB.fits'

DIFF_CCF_ORDER = 48
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

    # get files from ref path
    ref_path = os.path.join(paths[ref_name], obs_dir, filestr)
    reffiles = glob.glob(ref_path)

    filenames = []
    for reffile in reffiles:
        filenames.append(os.path.join(obs_dir, os.path.basename(reffile)))
    # ----------------------------------------------------------------------
    plt.close()
    fig, frames = plt.subplots(nrows=len(filenames), ncols=1,
                               figsize=(len(filenames)*4, 16))
    for it, filename in enumerate(filenames):

        frame = frames[it]

        data_ref = Table.read(os.path.join(paths[ref_name], filename), hdu=2)

        for name in paths:

            if name == ref_name:
                continue

            data_comp = Table.read(os.path.join(paths[name], filename), hdu=2)

            diff = 1000 * (data_comp['RV'] - data_ref['RV'])

            frame.plot(data_comp['Orders'], diff, label=f'{ref_name}-{name}',
                       marker='o', ls='None')

        frame.set(xlabel='Orders', ylabel='$\Delta$RV  m/s', title=filename)
        frame.legend(loc=0)

    plt.tight_layout()
    plt.show()
    plt.close()
    # ----------------------------------------------------------------------
    plt.close()
    fig, frames = plt.subplots(nrows=len(filenames), ncols=2,
                               figsize=(len(filenames)*4, 16))
    for it, filename in enumerate(filenames):
        frame1 = frames[it, 0]
        frame2 = frames[it, 1]

        data_ref = Table.read(os.path.join(paths[ref_name], filename), hdu=1)
        rv = data_ref['RV']

        for name in paths:

            if name == ref_name:
                continue

            data_comp = Table.read(os.path.join(paths[name], filename), hdu=1)

            o_ref = data_ref['ORDER{:02d}'.format(DIFF_CCF_ORDER)]
            o_comp = data_comp['ORDER{:02d}'.format(DIFF_CCF_ORDER)]
            diff = o_ref/o_comp
            frame1.plot(rv, o_ref, label=f'{ref_name}')
            frame1.plot(rv, o_comp, label=f'{name}')
            frame2.plot(rv, diff, label=f'{ref_name}-{name}')
        frame1.set(xlabel='RV km/s', ylabel=f'CCF {DIFF_CCF_ORDER}',
                  title=filename)
        frame2.set(xlabel='RV km/s', ylabel=f'ratio CCF {DIFF_CCF_ORDER}',
                  title=filename)
        frame1.legend(loc=0)
        frame2.legend(loc=0)
    plt.tight_layout()
    plt.show()
    plt.close()

# =============================================================================
# End of code
# =============================================================================
