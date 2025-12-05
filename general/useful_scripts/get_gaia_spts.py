#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-10-19 at 16:07

@author: cook
"""
from astroquery.gaia import Gaia
import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# Define variables
# =============================================================================
# The distance
DISTANCE = 10
# The SQL query
QUERY = """
SELECT a.*, b.* 
FROM gaiadr3.gaia_source AS a 
INNER JOIN gaiadr3.astrophysical_parameters as b ON a.source_id=b.source_id 
WHERE parallax > {parallax}
"""
# The unique spectral types
USPTS = ['O', 'B', 'A', 'F', 'G', 'K', 'M']
# the color for each spectral type
UCOLOURS = ['b', 'c', 'orange', 'purple', 'yellow', 'orange', 'red']
# The spectral type colomn
SPTCOL = "spectraltype_esphs"
# -----------------------------------------------------------------------------

# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # work out parallax in mas
    parallax = 1000 / DISTANCE
    # query gaia
    job = Gaia.launch_job_async(QUERY.format(parallax=parallax))
    # get the results
    table = job.get_results()

    # get the spectral type column data
    spt_data = table[SPTCOL]

    # store the number of each spectral type
    scounts = []
    # keep track of those found
    idmask = np.zeros(len(spt_data), dtype=bool)
    # count the number of each spectral type
    for uspt in USPTS:
        # find all rows of this spectral type
        mask = spt_data == uspt
        # add the count to the storage
        scounts.append(np.sum(mask))
        # keep track of those found
        idmask |= mask

    # identify those without spectral type above
    scounts.append(np.sum(~idmask))
    # add an unknown column
    USPTS.append('unknown')
    # add a colour
    UCOLOURS.append('k')
    # add an numerical x axis
    x = np.arange(len(scounts))
    # plot
    plt.close()
    fig, frame = plt.subplots()
    frame.bar(x, scounts, color=UCOLOURS)
    frame.set_xticks(x, USPTS)
    tkwargs = dict(dist=DISTANCE, N=len(table))
    title = 'Sample at {dist} pc [Total={N}]'.format(**tkwargs)
    frame.set(xlabel='Spectral Type', ylabel='Number of objects',
              title=title)
    plt.show()
    plt.close()


# =============================================================================
# End of code
# =============================================================================
