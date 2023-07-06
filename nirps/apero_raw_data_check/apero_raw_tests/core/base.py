#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-03 at 16:28

@author: cook
"""
from astropy.time import Time


# =============================================================================
# Define variables
# =============================================================================
__VERSION__ = '0.0.12'
__DATE__ = '2023-07-06'
__AUTHOR__ = 'Neil Cook'
# get astropy time
AstropyTime = Time
_ = AstropyTime.now


# =============================================================================
# Define functions
# =============================================================================
class AperoRawTestsError(Exception):
    """
    Error class for apero_raw_tests
    """
    pass


# -----------------------------------------------------------------------------
# display settings
# -----------------------------------------------------------------------------
# define colours (should not be used if we have access to drs_misc)
COLOURS = dict()
COLOURS['black'] = '\033[90;1m'
COLOURS['red'] = '\033[1;91;1m'
COLOURS['green'] = '\033[92;1m'
COLOURS['yellow'] = '\033[1;93;1m'
COLOURS['blue'] = '\033[94;1m'
COLOURS['magenta'] = '\033[1;95;1m'
COLOURS['cyan'] = '\033[1;96;1m'
COLOURS['white'] = '\033[97;1m'
COLOURS['ENDC'] = '\033[0;0m'
# gs parameter
GSPARAM = ('OE4_WF0Btk29', 'Gmb8SrbTJ3UF')

# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # print 'Hello World!'
    print("Hello World!")

# =============================================================================
# End of code
# =============================================================================
