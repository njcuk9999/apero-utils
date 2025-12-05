#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-06-16

@author: cook
"""
from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np


# =============================================================================
# Define variables
# =============================================================================
filename = '2022-06-11/NIRPS_2022-06-12T04_31_38_321_pp_s1d_v_A.fits'
# cmap = 'gist_ncar'
# cmap = 'gist_rainbow_r'
cmap = 'rainbow'




# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":

    table = Table.read(filename)

    mask = np.isfinite(table['wavelength']) & np.isfinite(table['flux'])

    x = table['wavelength'][mask]
    y = table['flux'][mask]
    z = x



    # Create a set of line segments so that we can color them individually
    # This creates the points as a N x 1 x 2 array so that we can stack points
    # together easily to get the segments. The segments array for line collection
    # needs to be (numlines) x (points per line) x 2 (for x and y)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    fig, frame = plt.subplots(1, 1)

    # Create a continuous norm to map from data points to colors
    norm = plt.Normalize(z.min(), z.max())
    lc = LineCollection(segments, cmap=cmap, norm=norm)
    # Set the values used for colormapping
    lc.set_array(z)
    lc.set_linewidth(0.5)
    line = frame.add_collection(lc)
    # fig.colorbar(line, ax=frame)

    # Use a boundary norm instead
    frame.set(xlabel='Wavelength [nm]', ylabel='Flux',
              xlim=[x.min(), x.max()], ylim=[0, y.max()])

    plt.show()

# =============================================================================
# End of code
# =============================================================================
