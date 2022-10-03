#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-10-03 at 15:28

@author: cook
"""
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

# =============================================================================
# Define variables
# =============================================================================
# Set the file name of the s
# FILENAME = '/nirps_raw/nirps/apero-data/drs-data/nirps_ha_202205/reduced/2022-06-18/NIRPS_2022-06-19T00_58_00_800_pp_s1d_v_A.fits'

FILE1 = '/nirps_raw/nirps/apero-data/drs-data/nirps_ha_202205/reduced/other/Template_s1d_PROXIMA_sc1d_v_file_A.fits'
FILE2 = '/nirps_raw/nirps/apero-data/drs-data/nirps_ha_202205/reduced/other/Template_s1d_GL699_sc1d_v_file_A.fits'
# -----------------------------------------------------------------------------

FILES = [FILE1, FILE2]


# =============================================================================
# Define functions
# =============================================================================


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    plt.close()
    fig, frames = plt.subplots(nrows=len(FILES), ncols=1, figsize=(30, 5*len(FILES)))


    for it, filename in enumerate(FILES):

        frame = frames[it]

        table = Table.read(filename, hdu=1)

        x = np.array(table['wavelength'])
        y = np.array(table['flux'])

        # mask values less than zero and greater than 1

        # mask = (x > 1370) & (x < 1405)

        # x[mask] = np.nan
        # y[mask] = np.nan

        # Create a set of line segments so that we can color them individually
        # This creates the points as a N x 1 x 2 array so that we can stack points
        # together easily to get the segments. The segments array for line collection
        # needs to be (numlines) x (points per line) x 2 (for x and y)
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)


        fig.patch.set_facecolor('k')
        frame.set_facecolor('k')
        # Create a continuous norm to map from data points to colors
        norm = plt.Normalize(np.nanmin(x), np.nanmax(x))
        lc = LineCollection(segments, cmap='autumn_r', norm=norm)
        # Set the values used for colormapping
        lc.set_array(x)
        lc.set_linewidth(0.75)
        lc.set_alpha(0.4)
        line = frame.add_collection(lc)

        frame.set(xlim=(np.nanmin(x), np.nanmax(x)),
                  ylim=(0, np.nanmax(y)))

        frame.tick_params(colors='white', which='both')

        frame.spines['bottom'].set_color('white')
        frame.spines['top'].set_color('white')
        frame.spines['right'].set_color('white')
        frame.spines['left'].set_color('white')
        frame.set_xticklabels([])
        frame.set_yticklabels([])
        frame.tick_params(axis='x', which='both', bottom=False, top=False)
        frame.tick_params(axis='y', which='both', left=False, right=False)

    plt.savefig('/data/spirou/cook/etienne_proxima.jpg', dpi=300)
    plt.close()

# =============================================================================
# End of code
# =============================================================================
