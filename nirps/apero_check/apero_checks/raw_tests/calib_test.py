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
from typing import Any, Dict, Tuple

from astropy.io import fits
from tqdm import tqdm

from apero_checks.core import misc
from apero_checks.core import io

# =============================================================================
# Define variables
# =============================================================================
# define any other constants here that you want moving to the parameters.py
#  file (these can be overwritten by the yaml file) and may change
#  depending on the profile used (i.e. NIRPS_HA or NIRPS_HE)


# =============================================================================
# Define functions
# =============================================================================
def test(params: Dict[str, Any], obsdir: str, log=False) -> Tuple[bool, str]:
    """
    Test with observation directory exists

    Passed = True

    :param params: dictionary of parameters (from yaml/cmd line/ parameters.py)
    :param obsdir: str, the observation directory (e.g. the night directory)
    :param log: bool, if True prints messages (all messages must be wrapped
                in a if log: statement)

    :return: bool, True if passed, False otherwise
    """

    # Required DPRTYPE groups. Each inner list is an OR condition.
    # For example, ['WAVE,UN1,FP', 'WAVE,UN2,FP'] means either is accepted.
    must_have_dprtypes = [
        ['DARK'],
        ['WAVE,FP,FP'],
        ['WAVE,UN1,FP', 'WAVE,UN2,FP'],
        ['WAVE,FP,UN1', 'WAVE,FP,UN2'],
        ['FLAT,DARK,LAMP'],
        ['FLAT,LAMP,DARK'],
        ['WAVE,UN1,UN1', 'WAVE,UN2,UN2']
    ]
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
        out_msg = ('Observation directory {0} does not exist - TEST FAILED')
        if log:
            print(out_msg)
        return False, out_msg
    # get the files for this observation directory
    files = glob.glob(os.path.join(obsdir_path, '*.fits'))
    # deal with no files
    if len(files) == 0:
        # pass a True if no file is found on that night
        out_msg = ('No files found for night {}'.format(obsdir))
        if log:
            print(out_msg)
        return False, out_msg
    # -------------------------------------------------------------------------
    # create table to store keywords
    dpr_counts = dict()
    # fill table looping through files
    for ifile in tqdm(range(len(files))):
        # get filename for this iteration
        filename = files[ifile]
        # get dpr type and obs name
        dpr_type = io.get_header_key(filename, header_keys['DPR_TYPE'],
                                     required=False, default=None)
        # count dpr type instances
        if dpr_type in dpr_counts:
            dpr_counts[dpr_type] += 1
        else:
            dpr_counts[dpr_type] = 1
    # -------------------------------------------------------------------------
    # set the failed messages
    failed_log = []
    passed_log = []
    dpr_type_qc_passed = True
    # -------------------------------------------------------------------------
    # now we test must have dpr types
    # -------------------------------------------------------------------------
    # could be updated to force a minimum number of calibrations
    for dpr_group in must_have_dprtypes:
        found = False
        matched_dpr = None
        for dpr_type in dpr_group:
            if dpr_type in dpr_counts:
                found = True
                matched_dpr = dpr_type
                break
        if not found:
            if len(dpr_group) == 1:
                fmsg = '\tNo "ESO DPR TYPE" = {} on that night'
                failed_log.append(fmsg.format(dpr_group[0]))
            else:
                dmsg = '\tNo "ESO DPR TYPE" in ({}) on that night'
                failed_log.append(dmsg.format(' OR '.join(dpr_group)))
            dpr_type_qc_passed = False
        else:
            if len(dpr_group) == 1:
                pmsg = '\tFound "ESO DPR TYPE" = {}'
                passed_log.append(pmsg.format(matched_dpr))
            else:
                pmsg = '\tFound "ESO DPR TYPE" = {} (from {})'
                pargs = [matched_dpr, ' OR '.join(dpr_group)]
                passed_log.append(pmsg.format(*pargs))
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
    outmsg = ('\n')
    outmsg += ('\n' + '*' * 50)
    outmsg += ('\nQC for night {}'.format(obsdir))
    outmsg += ('\nPassed QC:')
    outmsg += ('\n ' + passed_log)
    outmsg += ('\nFailed QC:')
    outmsg += ('\n' + failed_log)
    outmsg += ('\n' + '*' * 50)
    outmsg += ('\n')
    if log:
        print(outmsg)
    # -------------------------------------------------------------------------
    return passed_logical, outmsg


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
