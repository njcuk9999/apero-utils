#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-02-23 at 16:38

@author: cook
"""
from astropy.io import fits
import numpy as np

from apero.core.constants import constant_functions
from apero.core.instruments.spirou import default_keywords

# =============================================================================
# Define variables
# =============================================================================
# Keyword defintion
Keyword = constant_functions.Keyword
# Demo file
DEMO_FILE = '/spirou2/apero-data/common/minidata2/2020-08-31/2510309o.fits'
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
    # list of keys to check for in raw header
    raw_keys = []
    # get all attributes from default_keywords
    for key in dir(default_keywords):
        # get the attribute
        attr = getattr(default_keywords, key)
        # check if attribute is a apero keyword
        if isinstance(attr, Keyword):
            # add to list
            raw_keys.append(attr.key)
    # ----------------------------------------------------------------------
    # load the raw file header
    raw_header = fits.getheader(DEMO_FILE)
    # ----------------------------------------------------------------------
    # storage of demo key value and comment
    storage = dict()
    # find which raw keys are in demo file
    for key in raw_keys:
        if key in raw_header:
            storage[key] = (raw_header[key], raw_header.comments[key])

    # sort keys alphabetically
    sorted_keys = np.sort(list(storage.keys()))
    # print out the results in readable format
    print('Raw header keys found in demo file:')
    for key in sorted_keys:
        print('  {0:20s} = {1:20s} # {2}'.format(key, str(storage[key][0]),
                                                  storage[key][1]))




# =============================================================================
# End of code
# =============================================================================
