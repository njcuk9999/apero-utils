#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Blank example test do not modify this.

Copy this to create your test

Created on 2023-07-03 at 14:37

@author: cook
"""
import os
from typing import Any, Dict, Tuple
import pandas as pd
from astropy.time import Time


# =============================================================================
# Define variables
# =============================================================================
# define any other constants here that you want moving to the parameters.py
#  file (these can be overwritten by the yaml file) and may change
#  depending on the profile used (i.e. NIRPS_HA or NIRPS_HE)

# First test date (before this date there was no critical checks)
FIRST_TEST_DATE = Time('2025-08-25T00:00:00', format='fits')

# Define the types of critical tests
DESC_TYPES: Dict[str, str] = dict()
DESC_TYPES['CRITICAL'] = 'raw'
DESC_TYPES['CRITICAL_SCI'] = 'sci'


# =============================================================================
# Define functions
# =============================================================================
def critical_test(params: Dict[str, Any], obsdir: str, log=False
                  ) -> Tuple[bool, str]:
    return test(params, obsdir, 'CRITICAL', log=log)


def critical_sci_test(params: Dict[str, Any], obsdir: str, log=False
                      ) -> Tuple[bool, str]:
    return test(params, obsdir, 'CRITICAL_SCI', log=log)    


def test(params: Dict[str, Any], obsdir: str, tkind: str, 
         log=False) -> Tuple[bool, str]:
    """
    Blank test - this tests whether test was run (should always return True)
    All other tests should return True or False, and only print messages if
    log is True.

    Passed = True

    :param params: dictionary of parameters (from yaml/cmd line/ parameters.py)
    :param obsdir: str, the observation directory (e.g. the night directory)
    :param tkind: str, type of critical test (CRITICAL or CRITICAL_SCI)
    :param log: bool, if True prints messages (all messages must be wrapped
                in a if log: statement)

    :return: bool, True if passed, False otherwise
    """
    # -------------------------------------------------------------------------
    # get a machine readable date
    mr_obsdir = Time(obsdir + 'T00:00:00', format='fits')
    # -------------------------------------------------------------------------
    # Don't test if obsdir is before FIRST_TEST_DATE
    if mr_obsdir < FIRST_TEST_DATE:
        out_msg = ('\nCRITICAL TEST: Obsdir {0} is before first test date '
                   '{1}, automatically passing'.format(obsdir, FIRST_TEST_DATE))
        if log:
            print(out_msg)
        return True, out_msg
    # -------------------------------------------------------------------------
    # we don't use parameters here but all tests must take params as first
    #   and can use any argument from parameters
    csv_file = params['critical csv file']
    desc_file = params['critical desc file']
    # deal with override
    if csv_file == 'None' or desc_file == 'None':
        out_msg = 'No critial csv or desc file set. Automatically passing test.'

        print('"check.critical csv file = None" or '
              '"check.critical desc file = None" '
              '\n\t- No critical csv or desc file set, '
              'skipping test')

        return True, out_msg
    # -------------------------------------------------------------------------
    # files are only accessible on the NIRPS system - we get around this by
    #   checking ~/csv_file if csv_file does not exist (we assume testing
    #   user has drive mounted in home directory)
    if not os.path.exists(csv_file):
        alt_file = os.path.expanduser('~') + csv_file
        csv_file = alt_file
        desc_file = os.path.expanduser('~') + desc_file
        params['critical csv file'] = csv_file
        params['critical desc file'] = desc_file
    # -------------------------------------------------------------------------
    # Load the csv file using pandas
    try:
        df = pd.read_csv(csv_file, index_col=0)
    except Exception as _:
        out_msg = ('\nCRITICAL TEST: Could not read the critical checks '
                   'csv file: {0}'.format(csv_file))
        if log:
            print(out_msg)
        return False, out_msg
    # -------------------------------------------------------------------------
    # load the check descriptions
    try:
        df_desc = pd.read_csv(desc_file)
    except Exception as _:
        out_msg = ('\nCRITICAL TEST: Could not read the critical checks '
                   'description csv file: {0}'.format(desc_file))
        if log:
            print(out_msg)
        return False, out_msg

    # check that obsdir is in the index (if not return False)
    if obsdir not in df.index:
        out_msg = ('\nCRITICAL TEST: No entry for obsdir {0} in critical '
                   'checks csv file: {1}'.format(obsdir, csv_file))
        if log:
            print(out_msg)
        return False, out_msg
    # -------------------------------------------------------------------------
    # we have the obsdir so we can check the individual checks
    # -------------------------------------------------------------------------
    # get the row for this obsdir
    row = df.loc[obsdir]

    # loop around the columns (taken from df_desc) and check if any are False
    passed = True
    msgs = []
    for _, desc in df_desc.iterrows():
        # get the check name and description
        check_name = desc['name']
        check_desc = desc['description']
        check_type = desc['type']

        # if we have a type then check that it matches tkind
        # TODO: need to update the column name "type" once csv is updated
        # then we can remove the if 'type' part
        if 'type' in desc:
            check_type = desc['type'] 
            if check_type != DESC_TYPES[tkind]:
                continue
        # check if the check_name is in the row
        if check_name not in row:
            out_msg = ('\nCRITICAL TEST: Check name {0} not in critical '
                       'checks csv file: {1}'.format(check_name, csv_file))
            msgs.append(out_msg)
            passed = False
        # check if the row value is False
        elif not row[check_name]:
            out_msg = ('\nCRITICAL TEST: Check {0} failed - {1}'
                       .format(check_name, check_desc))
            msgs.append(out_msg)
            passed = False
        # otherwise the check passed
        else:
            out_msg = ('\nCRITICAL TEST: Check {0} passed - {1}'
                       .format(check_name, check_desc))
            msgs.append(out_msg)

    # combine output messages
    out_msgs = '\n'.join(msgs)
    # all print out messages must be wrapped in if log
    if log:
        print(out_msgs)
    # return whether passed and the messages
    return passed, out_msgs


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # provide a working example of a test
    _params = dict()
    # define the observation directory
    _obsdir = '2021-03-15'
    # run the test
    test(_params, _obsdir, log=True, tkind='CRITICAL')

# =============================================================================
# End of code
# =============================================================================
