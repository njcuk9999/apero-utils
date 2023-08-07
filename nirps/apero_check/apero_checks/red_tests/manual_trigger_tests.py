#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-08-07 at 11:44

@author: cook
"""

import os
from astropy.table import Table


# =============================================================================
# Define variables
# =============================================================================
# define messages
MANUAL_START = 'MANUAL_START'
MANUAL_END = 'MANUAL_END'
APERO_START = 'APERO_START'
APERO_ERR = 'APERO_ERR'
APERO_END = 'APERO_END'
LBL_START = 'LBL_START'
LBL_ERROR = 'LBL_ERR'
LBL_END = 'LBL_END'

MESSAGES = [MANUAL_START, MANUAL_END, APERO_START, APERO_ERR, APERO_END,
            LBL_START, LBL_ERROR, LBL_END]
# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def get_manual_log(profile_name, obs_dir):
    # construct log path
    homedir = os.path.expanduser('~')
    # construct log path
    logpath = os.path.join(homedir, '.apero', 'manual_trigger')
    # construct log file name
    logname = os.path.join(logpath, profile_name)
    # create a dictionary for tests
    tests = dict()
    # loop around tests and assume they are all False for now
    for test in MESSAGES:
        tests[test] = False


    # read csv file
    if not os.path.exists(logname):
        return tests

    # otherwise we read the log file
    log = Table.read(logname, format='csv')








# =============================================================================
# Define tests
# =============================================================================
def function1():
    return 0


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
