#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Blank example test do not modify this.

Copy this to create your test

Created on 2023-07-03 at 14:37

@author: cook
"""
import glob
import os
from typing import Any, Dict

from astropy.io import fits
from tqdm import tqdm

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

    # Here you can add more calibrations that are required (one day LFC maybe)
    must_have_dprtypes = ['DARK', 'WAVE,FP,FP', 'WAVE,UN1,FP', 'WAVE,FP,UN1',
                          'FLAT,DARK,LAMP', 'FLAT,LAMP,DARK', 'WAVE,UN1,UN1']
    # This tells use what the header keys are
    header_keys = dict()
    header_keys['DPR_TYPE'] = 'HIERARCH ESO DPR TYPE'
    header_keys['OBS_NAME'] = 'HIERARCH ESO OBS NAME'
    header_keys['MJD-OBS'] = 'MJD-OBS'
    # define parameters we use here
    raw_directory = params['raw dir']
    # -------------------------------------------------------------------------
    if log:
        msg = 'Analysing observation directory: {0}'
        margs = [obsdir]
        misc.log_msg(msg.format(*margs), level='')
    # -------------------------------------------------------------------------
    # get the observation directory
    obsdir_path = os.path.join(raw_directory, obsdir)
    # deal with no observation directory
    if not os.path.exists(obsdir_path):
        if log:
            print('Observation directory {0} does not exist - TEST FAILED')
        return False
    # get the files for this observation directory
    files = glob.glob(os.path.join(obsdir_path, '*.fits'))
    # deal with no files
    if len(files) == 0:
        if log:
            # pass a True if no file is found on that night
            print('No files found for night {}'.format(obsdir))
        return False
    # -------------------------------------------------------------------------
    # create table to store keywords
    dpr_counts = dict()
    # fill table looping through files
    for ifile in tqdm(range(len(files))):
        # get the header
        hdr = fits.getheader(files[ifile])
        # get dpr type and obs name
        dpr_type = str(hdr.get(header_keys['DPR_TYPE'], None))
        # count dpr type instances
        if dpr_type in dpr_counts:
            dpr_counts[dpr_type] += 1
        else:
            dpr_counts[dpr_type] = 1
        # close the header
        del hdr
    # -------------------------------------------------------------------------
    # set the failed messages
    failed_log = []
    passed_log = []
    dpr_type_qc_passed = True
    # -------------------------------------------------------------------------
    # now we test must have dpr types
    # -------------------------------------------------------------------------
    # could be updated to force a minimum number of calibrations
    for dpr_type in must_have_dprtypes:
        # deal with no values found for this night
        if dpr_type not in dpr_counts:
            fmsg = '\tNo "ESO DPR TYPE" = {} on that night'
            failed_log.append(fmsg.format(dpr_type))
            dpr_type_qc_passed = False
    # deal with passed log
    if dpr_type_qc_passed:
        passed_log.append('\tAll "must have" calibrations have been obtained '
                          'at least once')
    # -------------------------------------------------------------------------
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
    # -------------------------------------------------------------------------
    if log:
        print('\n')
        print('*' * 50)
        print('QC for night {}'.format(obsdir))
        print('Passed QC:')
        print(passed_log)
        print('Failed QC:')
        print(failed_log)
        print('*' * 50)
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
