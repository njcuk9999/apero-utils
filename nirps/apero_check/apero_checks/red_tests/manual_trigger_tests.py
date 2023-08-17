#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-08-07 at 11:44

@author: cook
"""
from typing import Any, Dict, Tuple

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
ARI_START = 'ARI_START'
ARI_END = 'ARI_END'


MESSAGES = [MANUAL_START, MANUAL_END, APERO_START, APERO_ERR, APERO_END,
            ARI_START, ARI_END]
# -----------------------------------------------------------------------------
# define the column names in the csv log file
COLUMNS = ['TIMESTAMP', 'PROFILE', 'STATUS', 'OBSDIRS', 'COMMENT']

# =============================================================================
# Define functions
# =============================================================================
def get_manual_log(params, obs_dir, test) -> Tuple[bool, str]:
    # construct log path
    homedir = os.path.expanduser('~')
    # construct log path
    logpath = os.path.join(homedir, '.apero', 'manual_trigger')
    # get profile name
    profile_name = os.path.basename(params['YAML']).replace('.yaml', '.log')
    # construct log file name
    logname = os.path.join(logpath, profile_name)
    # read csv file
    if not os.path.exists(logname):
        return False, 'Log name {0} does not exist'.format(logname)
    # otherwise we read the log file
    # no header line provided
    logtable = Table.read(logname, format='csv', data_start=0,
                          names=COLUMNS)
    # make all rows that have this message
    status_mask = logtable['STATUS'] == test
    # get cut down log table
    status_table = logtable[status_mask]
    # loop around rows in the logtable
    for row in range(len(status_table)):
        # get the list of observation directories for this row
        obs_dirs = status_table['OBSDIRS'][row].split('|')
        # if we have a * then skip
        if '*' in obs_dirs:
            continue
        if obs_dir in obs_dirs:
            return True, ''
    # if we get here return False
    return False, 'Not found in log {0}'.format(logname)


def test_switch(test_name: str, params: Dict[str, Any], obsdir: str,
                log=False) -> bool:
    # run test
    passed, reason = get_manual_log(params, obsdir, test=test_name)
    # log the failure if logging is required
    if log:
        if passed:
            print('{0} passed'.format(test_name))
        else:
            print('{0} failed: {1}'.format(test_name, reason))
    # return the passed condition
    return passed


# =============================================================================
# Define tests
# =============================================================================
def test_manual_trigger_start(params: Dict[str, Any], obsdir: str,
                              log=False) -> bool:
    # run the manual start test
    return test_switch(MANUAL_START, params, obsdir, log)


def test_manual_trigger_end(params: Dict[str, Any], obsdir: str,
                              log=False) -> bool:
    # run the manual start test
    return test_switch(MANUAL_END, params, obsdir, log)


def test_apero_start(params: Dict[str, Any], obsdir: str,
                              log=False) -> bool:
    # run the manual start test
    return test_switch(APERO_START, params, obsdir, log)


def test_apero_end(params: Dict[str, Any], obsdir: str,
                              log=False) -> bool:
    # run the manual start test
    return test_switch(APERO_END, params, obsdir, log)


def test_ari_start(params: Dict[str, Any], obsdir: str,
                              log=False) -> bool:
    # run the manual start test
    return test_switch(ARI_START, params, obsdir, log)


def test_ari_end(params: Dict[str, Any], obsdir: str,
                              log=False) -> bool:
    # run the manual start test
    return test_switch(ARI_END, params, obsdir, log)


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
