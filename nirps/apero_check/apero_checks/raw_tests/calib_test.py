#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Blank example test do not modify this.

Copy this to create your test

Created on 2023-07-03 at 14:37

@author: cook
"""
import copy
import glob
import os
from typing import Any, Dict

import numpy as np
from astropy.io import fits

from apero_checks.core import misc


# =============================================================================
# Define variables
# =============================================================================
# define any other constants here that you want moving to the parameters.py
#  file (these can be overwritten by the yaml file) and may change
#  depending on the profile used (i.e. NIRPS_HA or NIRPS_HE)


# =============================================================================
# Define functions
# =============================================================================
def test(params: Dict[str, Any], obsdir: str, log=False) -> bool:
    """
    Test with observation directory exists

    Passed = True

    :param params: dictionary of parameters (from yaml/cmd line/ parameters.py)
    :param obsdir: str, the observation directory (e.g. the night directory)
    :param log: bool, if True prints messages (all messages must be wrapped
                in a if log: statement)

    :return: bool, True if passed, False otherwise
    """
    # define parameters we use here
    raw_directory = params['raw dir']
    # -------------------------------------------------------------------------
    if log:
        msg = 'Analysing observation directory: {0}'
        margs = [obsdir]
        misc.log_msg(msg.format(*margs), level='')
    # -------------------------------------------------------------------------

    obsdir_path = os.path.join(raw_directory, obsdir)

    if not os.path.exists(obsdir_path):
        if log:
            print('Observation directory {0} does not exist - TEST FAILED')
        return False

    files = glob.glob(os.path.join(obsdir_path, '*.fits'))

    if len(files) == 0:
        if log:
            # pass a True if no file is found on that night
            print('No files found for night {}'.format(obsdir))
        return False

    # create table to store keywords
    tbl = dict()

    # keywords to extract and table name
    keys = [['DATE-OBS', 'DATE-OBS'],
            ['MJD-OBS', 'MJD-OBS'],
            ['HIERARCH ESO DPR TYPE', 'DPR_TYPE']]
    keys = np.array(keys)

    # define table columns
    hdr0 = fits.getheader(files[0])

    for i in range(len(keys[:, 0])):
        # loop on keys and which column to store them
        key = keys[i, 0]
        key2 = keys[i, 1]

        # if not present, we store a dummy value not to raise an error
        if key not in hdr0:
            val = 0.0
        else:
            val = hdr0[key]

        # pad strings to avoid error if there are long strings in later files
        if type(val) == str:
            val = ' ' * 100
        tbl[key2] = np.array([val] * len(files))

    # fill table looping through files
    for ifile in range(len(files)):
        hdr = fits.getheader(files[ifile])

        for i in range(len(keys[:, 0])):
            key = keys[i, 0]  # keyword in header
            key2 = keys[i, 1]  # keyword in table
            # only if key is present in header
            if key in hdr:
                tbl[key2][ifile] = copy.deepcopy(hdr[key])

    # change to True/False for strings
    for key in tbl.keys():
        if (tbl[key][0] == 'FALSE') or (tbl[key][0] == 'TRUE'):
            tbl[key] = tbl[key] == 'TRUE'

    # We add prints on the status of raw keywords:
    passed_log = []
    failed_log = []

    qc = True
    # Here you can add more calibrations that are required (one day LFC maybe)
    must_have_objects = ['DARK', 'WAVE,FP,FP', 'WAVE,UN1,FP', 'WAVE,FP,UN1', 'FLAT,DARK,LAMP',
                         'FLAT,LAMP,DARK', 'WAVE,UN1,UN1']

    for obj in must_have_objects:
        # could be updated to force a minimum number of calibrations
        qc = qc & np.any(tbl['DPR_TYPE'] == obj)

        if obj not in tbl['DPR_TYPE']:
            failed_log.append('\tNo "ESO DPR TYPE" = {} on that night'.format(obj))
    if qc:
        passed_log.append('\tAll "must have" calibrations have been obtained at least once')
    ####################################################################################

    # construct a string for printing output
    passed_log = '\n'.join(passed_log)
    if len(passed_log) == 0:
        passed_log = '\tNone'
    failed_log = '\n'.join(failed_log)

    if len(failed_log) == 0:
        failed_log = '\tNone'
        passed_logical = True
    else:
        passed_logical = False

    if log:
        print('\n')
        print('**********************************************************************')
        print('QC for night {}'.format(obsdir))
        print('Passed QC:')
        print(passed_log)
        print('Failed QC:')
        print(failed_log)
        print('**********************************************************************')
        print('\n')

    # -------------------------------------------------------------------------
    return passed_logical


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # provide a working example of a test
    _params = dict()
    # define the observation directory
    _obsdir = '2021-07-01'
    # run the test
    test(_params, _obsdir, log=True)

# =============================================================================
# End of code
# =============================================================================
